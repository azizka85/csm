#include <format>

#include <iostream>

#include <chrono>

#include <utils/bathymetry.h>
#include <utils/viscosity.h>
#include <utils/data.h>

#include "kernel.h"
#include "solver.h"

using namespace std::chrono;

using namespace Parallel::WindInducedCurrents::Davies85;

Solver::Solver(
	PhysicalParams phys,
	double g, double kb,
	double dx, double dy,
	double endTime, double outputTimeStep, string outDir,
	shared_ptr<Calc::Time::IAdjustTimeStep> timeStep,
	shared_ptr<Generators::Area::IGenerator> areaGenerator,
	shared_ptr<Generators::DZ::IGenerator> dzGenerator,
	shared_ptr<Generators::Bathymetry::IGenerator> hGenerator,
	shared_ptr<Generators::Wind::IGenerator> qGenerator,
	shared_ptr<Generators::Viscosity::IGenerator> nuGenerator
) {
	this->phys = phys;

	this->g = g;
	this->kb = kb;

	this->dx = dx;
	this->dy = dy;

	this->endTime = endTime;
	this->outputTimeStep = outputTimeStep;
	this->outDir = outDir;

	this->timeStep = timeStep;
	this->areaGenerator = areaGenerator;
	this->dzGenerator = dzGenerator;
	this->hGenerator = hGenerator;
	this->qGenerator = qGenerator;
	this->nuGenerator = nuGenerator;
}

Directories Solver::createDirectory() {
	auto dirPath = path(
		format("{}/f={}, g={}, rho={}, kb={}", outDir, phys.f, g, phys.rho, kb)
	);

	create_directories(dirPath);

	dirPath = hGenerator->addDirectory(dirPath);

	create_directories(dirPath);

	dirPath = nuGenerator->addDirectory(dirPath);

	create_directories(dirPath);

	dirPath = qGenerator->addDirectory(dirPath);

	create_directories(dirPath);

	dirPath = areaGenerator->addDirectory(dirPath);

	create_directories(dirPath);

	dirPath = dirPath / format("dx={}, dy={}", dx, dy);

	create_directories(dirPath);

	dirPath = dzGenerator->addDirectory(dirPath);

	create_directories(dirPath);

	auto surfacePath = dirPath / path("surface");

	create_directories(surfacePath);

	auto volumePath = dirPath / path("volume");

	create_directories(volumePath);

	auto viscosityPath = dirPath / path("viscosity");

	create_directories(viscosityPath);

	return Directories{
		.root = dirPath,
		.surface = surfacePath,
		.volume = volumePath,
		.viscosity = viscosityPath
	};
}

void Solver::setInitialCondition(
	int nx, int ny, int nz,
	vector<double>& uf,
	vector<double>& vf,
	vector<double>& ua,
	vector<double>& va,
	vector<double>& z
) {
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			int p = j + i * ny;

			for (int k = 0; k < nz; k++) {
				int id = k + p * nz;

				uf[id] = 0;
				vf[id] = 0;
			}

			ua[p] = 0;
			va[p] = 0;
			z[p] = 0;
		}
	}
}

void Solver::writeData(
	const TimeState& time,
	const VolumeState& volume,
	Current& current,
	Wind& currentWindData,
	vector<double>& z,
	const Directories& dirs,
	vector<double>& uf,
	vector<double>& ua,
	vector<double>& u,
	vector<double>& vf,
	vector<double>& va,
	vector<double>& v,
	vector<double>& w,
	vector<double>& nu,
	Writers::Surface::Writer& surfaceWriter,
	Writers::Volume::Writer& volumeElemWriter,
	Writers::Volume::Writer& volumeNodeWriter
) {
	surfaceWriter.write(
		TimeState{ .m = time.m, .t = time.t },
		SurfaceState{ .nx = volume.nx, .ny = volume.ny, .dx = volume.dx, .dy = volume.dy },
		vector<Column>{
		Column{ .name = "va", .precision = 3 },
			Column{ .name = "w10", .precision = 3 },
			Column{ .name = "q", .precision = 3 },
			Column{ .name = "z", .precision = 7 }
	},
		vector<Writers::Surface::DataKind>{
		SurfaceVector{
			.vx = current.ua,
			.vy = current.va
		},
			SurfaceVector{
				.vx = currentWindData.speed.u10,
				.vy = currentWindData.speed.v10,
		},
		SurfaceVector{
			.vx = currentWindData.stress.qx,
			.vy = currentWindData.stress.qy
		},
			z
	},
		dirs.surface
	);

	Utils::Data::volumePlusSurface(volume.nx, volume.ny, volume.nz, uf, ua, u);
	Utils::Data::volumePlusSurface(volume.nx, volume.ny, volume.nz, vf, va, v);

	volumeElemWriter.write(
		TimeState{ .m = time.m, .t = time.t },
		VolumeState{ .nx = volume.nx, .ny = volume.ny, .nz = volume.nz, .dx = volume.dx, .dy = volume.dy, .dz = volume.dz, .h = volume.h },
		vector<Column>{
		Column{ .name = "v", .precision = 3 }
	},
		vector<Writers::Volume::DataKind>{
		VolumeVector{ .vx = u, .vy = v, .vz = w }
	},
		dirs.volume
	);

	volumeNodeWriter.write(
		TimeState{ .m = time.m, .t = time.t },
		VolumeState{ .nx = volume.nx, .ny = volume.ny, .nz = volume.nz + 1, .dx = volume.dx, .dy = volume.dy, .dz = volume.dz, .h = volume.h },
		vector<Column>{
		Column{ .name = "nu", .precision = 3 }
	},
		vector<Writers::Volume::DataKind>{
		nu
	},
		dirs.viscosity
	);
}

void Solver::solve() {
	auto dirs = createDirectory();

	auto geo = areaGenerator->generate();

	size_t nx = static_cast<int>(
		ceil(geo.l / dx)
		) + 1;

	size_t ny = static_cast<int>(
		ceil(geo.w / dy)
		) + 1;

	auto dz = dzGenerator->generate();

	size_t nz = dz.size();

	auto h = hGenerator->generate(nx, ny);

	auto [hMin, hMax] = Utils::Bathymetry::minMaxH(h, nx, ny);

	double dzMin = *min_element(dz.begin(), dz.end());

	const double dtMax = min(dx, dy) / sqrt(2 * g * hMax) / 1.5;

	vector<double> ua(nx * ny);
	vector<double> va(nx * ny);

	vector<double> z(nx * ny);

	vector<double> uf(nx * ny * nz);
	vector<double> vf(nx * ny * nz);

	vector<double> u(nx * ny * nz);
	vector<double> v(nx * ny * nz);
	vector<double> w(nx * ny * nz);

	setInitialCondition(
		nx, ny, nz,
		uf, vf,
		ua, va,
		z
	);

	auto windData = qGenerator->generate(nx, ny);

	size_t currentWindIndex = 0;

	vector<double> u10 = windData[currentWindIndex].speed.u10;
	vector<double> v10 = windData[currentWindIndex].speed.v10;

	vector<double> qx = windData[currentWindIndex].stress.qx;
	vector<double> qy = windData[currentWindIndex].stress.qy;

	auto nu = nuGenerator->generate(
		VolumeState{
			.nx = nx, .ny = ny, .nz = nz + 1,
			.dx = dx, .dy = dy,
			.dz = dz, .h = h
		},
		WindCurrent{
			.speed = WindSpeed {
				.u10 = u10,
				.v10 = v10
			},
			.stress = WindStress {
				.qx = qx,
				.qy = qy
			},
			.current = Current {
				.ua = ua,
				.va = va
			}
		}
	);

	double nuMax = Utils::Viscosity::maxNU(nu, nx, ny, nz + 1);

	double dt = min(dtMax, hMin * dzMin * hMin * dzMin / 2 / nuMax);
	double dtp = dt;

	bool updateDT = false;

	double t = 0;
	double tn = outputTimeStep;

	size_t n = 1;
	size_t m = 0;

	vector<vector<double>> statistics;

	Writers::Statistics::Writer statisticsWriter;

	Writers::Surface::Writer heightWriter("heights.vtk", Writers::Utils::defaultWriteSurfaceDataHeader);
	Writers::Surface::Writer surfaceWriter("", Writers::Utils::defaultWriteSurfaceDataHeader);

	Writers::Volume::Writer volumeElemWriter("", Writers::Utils::writeVolumeElemHeader);
	Writers::Volume::Writer volumeNodeWriter("", Writers::Utils::writeVolumeNodeHeader);

	heightWriter.write(
		TimeState{ .m = m, .t = t },
		SurfaceState{ .nx = nx, .ny = ny, .dx = dx, .dy = dy }, {
			Column {.name = "h", .precision = 3 }
		}, {
			h
		},
		dirs.root
	);

	VolumeState volume = { .nx = nx, .ny = ny, .nz = nz, .dx = dx, .dy = dy, .dz = dz, .h = h };

	Current current = { .ua = ua, .va = va };

	writeData(
		TimeState{ .m = m, .t = t },
		volume,
		current,
		windData[currentWindIndex],
		z,
		dirs,
		uf, ua, u,
		vf, va, v,
		w, nu,
		surfaceWriter,
		volumeElemWriter,
		volumeNodeWriter
	);

	Kernel kernel(
		nx, ny, 
		dx, dy, 
		phys.rho, phys.f, kb, g, 
		dz, h, 
		ua, va, z, 
		uf, vf, 
		u10, v10,
		qx, qy,
		nu
	);

	long long calcTime = 0;

	auto start = high_resolution_clock::now();

	while (t <= endTime) {
		kernel.calcUA(dt);
		kernel.swapUA();

		kernel.calcVA(dt);
		kernel.swapVA();

		kernel.calcRHS(dt);
		kernel.createTridiagonalMatrix(dt);
		kernel.calcUVF();

		kernel.calcZ(dt);

		t += dt;
		n += 1;

		if (currentWindIndex < windData.size() - 1 && t >= windData[currentWindIndex + 1].time) {
			currentWindIndex += 1;

			u10 = windData[currentWindIndex].speed.u10;
			v10 = windData[currentWindIndex].speed.v10;

			qx = windData[currentWindIndex].stress.qx;
			qy = windData[currentWindIndex].stress.qy;

			kernel.updateQUV10Device();

			updateDT = true;
		}

		nu = nuGenerator->generate(
			VolumeState{
				.nx = nx, .ny = ny, .nz = nz + 1,
				.dx = dx, .dy = dy,
				.dz = dz, .h = h
			},
			WindCurrent{
				.speed = WindSpeed {
					.u10 = u10,
					.v10 = v10
				},
				.stress = WindStress {
					.qx = qx,
					.qy = qy
				},
				.current = Current {
					.ua = ua,
					.va = va
				}
			}
		);

		nuMax = Utils::Viscosity::maxNU(nu, nx, ny, nz + 1);

		double dt1 = min(dtMax, hMin * dzMin * hMin * dzMin / 2 / nuMax);

		if (dt1 < dtp || updateDT) {
			dtp = dt1;

			dt = timeStep->calculate(Time::State{
				.t = t,
				.dt = dt1,
				.dtMax = dtMax,
				.mult = false
			});

			updateDT = false;

			cout << format("Update dt = {}", dt) << endl;
		}
		else {
			dt = timeStep->calculate(Time::State{
				.t = t,
				.dt = dt,
				.dtMax = dtMax,
				.mult = true
			});
		}

		if (t >= tn) {
			auto end = high_resolution_clock::now();

			auto duration = duration_cast<milliseconds>(end - start).count();

			calcTime += duration;

			kernel.updateUVA();
			kernel.updateUVF();
			kernel.updateZ();
			kernel.updateNU();

			writeData(
				TimeState{ .m = m, .t = t },
				volume,
				current,
				windData[currentWindIndex],
				z,
				dirs,
				uf, ua, u,
				vf, va, v,
				w, nu,
				surfaceWriter,
				volumeElemWriter,
				volumeNodeWriter
			);

			auto umd = kernel.maxAbsUFDifference();
			auto vmd = kernel.maxAbsVFDifference();

			auto zmd = kernel.maxAbsZDifference();

			auto maxRes = kernel.maxAbsResidual();

			cout << format(
				"Write data in file t={:.3f}, nu={:.5f}, maxAbsRes={}, convergence of u={:.5f}, v={:.5f}, z={:.9f} with dt={:.5}, calc time={}",
				t, nuMax, maxRes, umd, vmd, zmd, dt, calcTime / 1000
			) << endl;

			statistics.push_back({ (double)n, tn, (double)(calcTime / 1000), maxRes, umd, vmd, zmd });

			kernel.updateUVFTDevice();
			kernel.updateZTDevice();

			m += 1;
			tn = t + outputTimeStep;

			start = high_resolution_clock::now();
		}
	}

	statisticsWriter.write(
		vector<Writers::Statistics::Column>{
		Writers::Statistics::Column{ .name = "n", .precision = 0 },
			Writers::Statistics::Column{ .name = "t", .precision = 3 },
			Writers::Statistics::Column{ .name = "calc_time", .precision = 0 },
			Writers::Statistics::Column{ .name = "max_res", .precision = 18 },
			Writers::Statistics::Column{ .name = "umd", .precision = 5 },
			Writers::Statistics::Column{ .name = "vmd", .precision = 5 },
			Writers::Statistics::Column{ .name = "zmd", .precision = 9 }
	},
		statistics,
		dirs.root
	);

	cout << "Total number of iterations: " << n;
}
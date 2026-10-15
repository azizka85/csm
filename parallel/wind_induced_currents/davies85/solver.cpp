#include <format>

#include <iostream>

#include <chrono>

#include <utils/bathymetry.h>
#include <utils/viscosity.h>
#include <utils/data.h>

#include "solver.h"

using namespace std::chrono;

using namespace Parallel::WindInducedCurrents::Davies85;

Solver::Solver(
	Data& data, 
	bool calcNU, bool interpolateWind,
	double endTime, double outputTimeStep, string outDir,
	shared_ptr<Calc::Time::IAdjustTimeStep> timeStep,
	shared_ptr<IKernel> kernel
): data(data) {
	this->calcNU = calcNU;
	this->interpolateWind = interpolateWind;

	this->endTime = endTime;
	this->outputTimeStep = outputTimeStep;
	this->outDir = outDir;

	this->timeStep = timeStep;	
	this->kernel = kernel;
}

Directories Solver::createDirectory() {
	auto dirPath = path(
		format("{}/f={}, g={}, rho={}, kb={}", outDir, data.f, data.g, data.rho, data.kb)
	);

	dirPath = kernel->addDirectory(dirPath);
	dirPath = data.addDirectory(dirPath);

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
	const WindCurrent& currentWindData,
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
	size_t nx = data.symmetrise ? volume.nx / 2 : volume.nx;
	size_t ny = data.symmetrise ? volume.ny / 2 : volume.ny;

	surfaceWriter.write(
		TimeState{ .m = time.m, .t = time.t },
		SurfaceState{ .nx = nx, .ny = ny, .bathymetry = volume.bathymetry },
		vector<Column>{
			Column{ .name = "va", .precision = 3 },
			Column{ .name = "w10", .precision = 3 },
			Column{ .name = "q", .precision = 3 },
			Column{ .name = "z", .precision = 7 }
		},
		vector<Writers::Surface::DataKind>{
			SurfaceVector{
				.vx = currentWindData.current.ua,
				.vy = currentWindData.current.va
			},
			SurfaceVector{
				.vx = currentWindData.speed.u10,
				.vy = currentWindData.speed.v10,
			},
			SurfaceVector{
				.vx = currentWindData.stress.qx,
				.vy = currentWindData.stress.qy
			},
			Scalar{
				.v = z
			}
		},
		dirs.surface
	);

	Utils::Data::volumePlusSurface(volume.nx, volume.ny, volume.nz, uf, ua, u);
	Utils::Data::volumePlusSurface(volume.nx, volume.ny, volume.nz, vf, va, v);

	volumeElemWriter.write(
		TimeState{ .m = time.m, .t = time.t },
		VolumeState{ .nx = nx, .ny = ny, .nz = volume.nz, .dz = volume.dz, .bathymetry = volume.bathymetry },
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
		VolumeState{ .nx = nx, .ny = ny, .nz = volume.nz + 1, .dz = volume.dz, .bathymetry = volume.bathymetry },
		vector<Column>{
			Column{ .name = "nu", .precision = 3 }
		},
		vector<Writers::Volume::DataKind>{
			Scalar{
				.v = nu
			}
		},
		dirs.viscosity
	);
}

void Solver::solve() {
	auto dirs = createDirectory();	

	auto [hMin, hMax] = Utils::Bathymetry::minMaxH(data.bathymetry.h, data.nx, data.ny);

	double dzMin = *min_element(data.dz.begin(), data.dz.end());

	const double dtMax = min(data.dx, data.dy) / sqrt(2 * data.g * hMax) / 1.5;	

	vector<double> u(data.nx * data.ny * data.nz);
	vector<double> v(data.nx * data.ny * data.nz);
	vector<double> w(data.nx * data.ny * data.nz);

	setInitialCondition(
		data.nx, data.ny, data.nz,
		data.uf, data.vf,
		data.ua, data.va,
		data.z
	);

	WindSpeed windSpeed{
		.u10 = data.windData[data.currentWindIndex].speed.u10,
		.v10 = data.windData[data.currentWindIndex].speed.v10
	};

	WindStress windStress{
		.qx = data.windData[data.currentWindIndex].stress.qx,
		.qy = data.windData[data.currentWindIndex].stress.qy
	};

	double nuMax = Utils::Viscosity::maxNU(data.nu, data.nx, data.ny, data.nz + 1);

	double dt = min(dtMax, hMin * dzMin * hMin * dzMin / 2 / nuMax);
	double dtp = dt;

	bool updateDT = false;

	auto zp = data.z;

	auto up = data.uf;
	auto vp = data.vf;

	double t = 0;
	double tn = outputTimeStep;

	double t0 = data.windData[data.currentWindIndex].time;
	double t1 = data.windData[data.currentWindIndex].time;

	if (interpolateWind && data.currentWindIndex + 1 < data.windData.size()) {
		t1 = data.windData[data.currentWindIndex + 1].time;		
	}

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
		SurfaceState{ .nx = data.symmetrise ? data.nx / 2 : data.nx, .ny = data.symmetrise ? data.ny / 2 : data.ny, .bathymetry = data.bathymetry }, 
		{
			Column {.name = "h", .precision = 3 }
		}, 
		{
			Scalar {
				.v = data.bathymetry.h
			}
		},
		dirs.root
	);

	VolumeState volume = { .nx = data.nx, .ny = data.ny, .nz = data.nz, .dz = data.dz, .bathymetry = data.bathymetry };

	Current current{
		.ua = data.ua,
		.va = data.va
	};

	writeData(
		TimeState{ .m = m, .t = t },
		volume,
		WindCurrent{
			.speed = windSpeed,
			.stress = windStress,
			.current = current
		},
		data.z,
		dirs,
		data.uf, data.ua, u,
		data.vf, data.va, v,
		w, data.nu,
		surfaceWriter,
		volumeElemWriter,
		volumeNodeWriter
	);

	long long calcTime = 0;

	auto start = high_resolution_clock::now();

	m += 1;

	while (t <= endTime) {
		kernel->calcUA(dt);
		kernel->swapUA();

		kernel->calcVA(dt);
		kernel->swapVA();

		kernel->calcRHS(dt);
		kernel->createTridiagonalMatrix(dt);
		kernel->calcUVF();

		kernel->calcZ(dt);

		t += dt;
		n += 1;

		if (data.currentWindIndex + 1 < data.windData.size() && t >= data.windData[data.currentWindIndex + 1].time) {
			data.currentWindIndex += 1;

			kernel->updateQUV10Device(interpolateWind);

			if (interpolateWind && data.currentWindIndex + 1 < data.windData.size()) {				
				t0 = data.windData[data.currentWindIndex].time;
				t1 = data.windData[data.currentWindIndex + 1].time;
			}

			updateDT = true;
		}

		if (interpolateWind) {
			double alpha = (t1 == t0) ? 0 : (t - t0) / (t1 - t0);

			kernel->interpolateQUV10(alpha);
		}

		if (calcNU) {
			kernel->calcNU();
			nuMax = kernel->calcNUMax();
		}

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

			std::cout << format("Update dt = {}", dt) << endl;
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

			kernel->calcResidual();

			kernel->updateUVA();
			kernel->updateUVF();
			kernel->updateUVD();
			kernel->updateZ();

			if (calcNU) {
				kernel->updateNU();
			}

			if (interpolateWind) {
				kernel->updateQUV10(windSpeed.u10, windSpeed.v10, windStress.qx, windStress.qy);
			}

			writeData(
				TimeState{ .m = m, .t = t },
				volume,
				WindCurrent{
					.speed = windSpeed,
					.stress = windStress,
					.current = current
				},
				data.z,
				dirs,
				data.uf, data.ua, u,
				data.vf, data.va, v,
				w, data.nu,
				surfaceWriter,
				volumeElemWriter,
				volumeNodeWriter
			);

			auto umd = Utils::Data::maxAbsDifference(data.nx, data.ny, data.nz, up, data.uf);
			auto vmd = Utils::Data::maxAbsDifference(data.nx, data.ny, data.nz, vp, data.vf);
			auto zmd = Utils::Data::maxAbsDifference(data.nx, data.ny, zp, data.z);

			auto maxRes = max(
				Utils::Data::maxAbsElem(data.nx, data.ny, data.nz, data.ud),
				Utils::Data::maxAbsElem(data.nx, data.ny, data.nz, data.vd)
			);

			cout << format(
				"Write data in file t={:.3f}, nu={:.5f}, maxAbsRes={}, convergence of u={:.5f}, v={:.5f}, z={:.19f} with dt={:.5}, calc time={}",
				t, nuMax, maxRes, umd, vmd, zmd, dt, calcTime / 1000
			) << endl;

			statistics.push_back({ (double)n, tn, (double)(calcTime / 1000), maxRes, umd, vmd, zmd });

			up = data.uf;
			vp = data.vf;
			
			zp = data.z;

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
			Writers::Statistics::Column{ .name = "zmd", .precision = 19 }
	},
		statistics,
		dirs.root
	);

	cout << "Total number of iterations: " << n;
}
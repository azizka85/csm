#include <format>

#include <iostream>

#include <chrono>

#include <slae/direct/tridiagonal.h>

#include <utils/bathymetry.h>
#include <utils/viscosity.h>
#include <utils/data.h>

#include "solver.h"

using namespace std::chrono;

using namespace SLAE::Direct;

using namespace WindInducedCurrents::Davies85;

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

void Solver::calcUA(
	int nx, int ny, int nz, 
	PhysicalParams phys, 
	double kb, double g, 
	double dx, double dt, 
	vector<double>& h, 
	vector<double>& qx, 
	vector<double>& ua, 
	vector<double>& va, 
	vector<double>& z, 
	vector<double>& uf, 
	vector<double>& u1a
) {
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			size_t k = j + i * ny;
			size_t idx = nz - 1 + nz * k;

			if (h[k] > 0) {
				u1a[k] = ua[k] + phys.f * va[k] * dt - kb * dt * (uf[idx] + ua[k]) / h[k] + qx[k] * dt / phys.rho / h[k];

				if (i > 0 && i < nx - 1) {
					int kr = j + (i + 1) * ny;
					int kl = j + (i - 1) * ny;


					u1a[k] -= g * dt * (z[kr] - z[kl]) / 2 / dx;
				}
			}
		}
	}
}

void Solver::calcVA(
	int nx, int ny, int nz, 
	PhysicalParams phys, 
	double kb, double g, 
	double dy, double dt, 
	vector<double>& h, 
	vector<double>& qy, 
	vector<double>& ua, 
	vector<double>& va, 
	vector<double>& z, 
	vector<double>& vf, 
	vector<double>& v1a
) {
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			size_t k = j + i * ny;
			size_t idx = nz - 1 + nz * k;

			if (h[k] > 0) {
				v1a[k] = va[k] - phys.f * ua[k] * dt - kb * dt * (vf[idx] + va[k]) / h[k] + qy[k] * dt / phys.rho / h[k];

				if (j > 0 && j < ny - 1) {
					int ku = j + 1 + i * ny;
					int kd = j - 1 + i * ny;

					v1a[k] -= g * dt * (z[ku] - z[kd]) / 2 / dy;
				}
			}
		}
	}
}

void Solver::calcRHS(
	int nx, int ny, int nz, 
	PhysicalParams phys, double kb, 
	double dx, double dy, double dt, 
	vector<double>& dz, 
	vector<double>& h, 
	vector<double>& nu, 
	vector<double>& qx, 
	vector<double>& qy, 
	vector<double>& ua, 
	vector<double>& u1a, 
	vector<double>& va, 
	vector<double>& v1a, 
	vector<double>& uf, 
	vector<double>& vf, 
	vector<double>& ud, 
	vector<double>& vd
) {
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			size_t p = j + i * ny;

			if (h[p] > 0) {
				for (size_t k = 0; k < nz; k++) {
					size_t id = k + p * nz;

					int idu = k + 1 + p * nz;
					int idd = k - 1 + p * nz;

					int idb = nz - 1 + p * nz;

					if (k == 0) {
						double um1 = uf[id] + h[p] * qx[p] * dz[k] / phys.rho / nu[id];
						double vm1 = vf[id] + h[p] * qy[p] * dz[k] / phys.rho / nu[id];

						ud[id] = uf[id] + phys.f * vf[id] * dt + dt * (
							nu[idu] * (uf[idu] - uf[id]) / (dz[k + 1] + dz[k]) -
							nu[id] * (uf[id] - um1) / 2 / dz[k]
						) / h[p] / h[p] / dz[k] + kb * (uf[idb] + ua[p]) * dt / h[p] - qx[p] * dt / phys.rho / h[p]
						+ dt * qx[p] / phys.rho / h[p] / 2 / dz[k];

						vd[id] = vf[id] - phys.f * uf[id] * dt + dt * (
							nu[idu] * (vf[idu] - vf[id]) / (dz[k + 1] + dz[k]) -
							nu[id] * (vf[id] - vm1) / 2 / dz[k]
						) / h[p] / h[p] / dz[k] + kb * (vf[idb] + va[p]) * dt / h[p] - qy[p] * dt / phys.rho / h[p]
						+ dt * qy[p] / phys.rho / h[p] / 2 / dz[k];
					}
					else if (k == nz - 1) {
						double un = (
							uf[id] - kb * h[p] * dz[k] * uf[id] / 2 / nu[idu]
							- kb * h[p] * dz[k] * ua[p] / nu[idu]
						) / (1 + kb * h[p] * dz[k] / 2 / nu[idu]);

						double vn = (
							vf[id] - kb * h[p] * dz[k] * vf[id] / 2 / nu[idu]
							- kb * h[p] * dz[k] * va[p] / nu[idu]
						) / (1 + kb * h[p] * dz[k] / 2 / nu[idu]);

						ud[id] = uf[id] + phys.f * vf[id] * dt + dt * (
							nu[idu] * (un - uf[id]) / 2 / dz[k] -
							nu[id] * (uf[id] - uf[idd]) / (dz[k] + dz[k - 1])
						) / h[p] / h[p] / dz[k] + kb * (uf[idb] + ua[p]) * dt / h[p] - qx[p] * dt / phys.rho / h[p]
						- dt * kb * ua[p] / h[p] / 2 / dz[k] / (1 + kb * h[p] * dz[k] / 2 / nu[idu]);

						vd[id] = vf[id] - phys.f * uf[id] * dt + dt * (
							nu[idu] * (vn - vf[id]) / 2 / dz[k] -
							nu[id] * (vf[id] - vf[idd]) / (dz[k] + dz[k - 1])
						) / h[p] / h[p] / dz[k] + kb * (vf[idb] + va[p]) * dt / h[p] - qy[p] * dt / phys.rho / h[p]
						- dt * kb * va[p] / h[p] / 2 / dz[k] / (1 + kb * h[p] * dz[k] / 2 / nu[idu]);
					}
					else {
						double u = uf[idb] + ua[p];
						double v = vf[idb] + va[p];

						ud[id] = uf[id] + phys.f * vf[id] * dt + dt * (
							nu[idu] * (uf[idu] - uf[id]) / (dz[k + 1] + dz[k]) -
							nu[id] * (uf[id] - uf[idd]) / (dz[k] + dz[k - 1])
						) / h[p] / h[p] / dz[k] + kb * u * dt / h[p] - qx[p] * dt / phys.rho / h[p];

						vd[id] = vf[id] - phys.f * uf[id] * dt + dt * (
							nu[idu] * (vf[idu] - vf[id]) / (dz[k + 1] + dz[k]) -
							nu[id] * (vf[id] - vf[idd]) / (dz[k] + dz[k - 1])
						) / h[p] / h[p] / dz[k] + kb * v * dt / h[p] - qy[p] * dt / phys.rho / h[p];
					}
				}
			}					
		}
	}
}

void Solver::createTridiagonalMatrix(
	double kb, double dt, 
	int nx, int ny, int nz, 
	vector<double>& dz, 
	vector<double>& h, 
	vector<double>& nu, 
	vector<double>& al, 
	vector<double>& ac, 
	vector<double>& ar
) {
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			size_t p = j + i * ny;

			if (h[p] > 0) {
				for (size_t k = 0; k < nz; k++) {
					size_t id = k + p * nz;
					size_t idu = k + 1 + p * nz;

					if (k == 0) {
						ac[id] = 1 + dt * nu[idu] / h[p] / h[p] / dz[k] / (dz[k + 1] + dz[k]);
						ar[id] = -dt * nu[idu] / h[p] / h[p] / dz[k] / (dz[k + 1] + dz[k]);
					}
					else if (k == nz - 1) {
						al[id] = -dt * nu[id] / h[p] / h[p] / dz[k] / (dz[k] + dz[k - 1]);
						ac[id] = 1 + dt * (
							kb * h[p] * dz[k] / (1 + kb * h[p] * dz[k] / 2 / nu[idu]) / 2 / dz[k] +
							nu[id] / (dz[k] + dz[k - 1])
						) / h[p] / h[p] / dz[k];
					}
					else {
						al[id] = -dt * nu[id] / h[p] / h[p] / dz[k] / (dz[k] + dz[k - 1]);
						ac[id] = 1 + dt * (
							nu[idu] / (dz[k + 1] + dz[k]) +
							nu[id] / (dz[k] + dz[k - 1])
						) / h[p] / h[p] / dz[k];
						ar[id] = -dt * nu[idu] / h[p] / h[p] / dz[k] / (dz[k + 1] + dz[k]);
					}
				}
			}
		}
	}
}

void Solver::calcUVF(
	int nx, int ny, int nz, 
	vector<double>& h, 
	vector<double>& al, 
	vector<double>& ac, 
	vector<double>& ar, 
	vector<double>& uf, 
	vector<double>& vf, 
	vector<double>& ud, 
	vector<double>& vd
) {
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			size_t k = j + i * ny;			

			if (h[k] > 0) {
				Tridiagonal::Matrix matrix(
					Tridiagonal::MatrixParams{
						.l = span<double>(al.data() + k * nz, nz),
						.c = span<double>(ac.data() + k * nz, nz),
						.r = span<double>(ar.data() + k * nz, nz)
					}
				);

				matrix.solve(
					Tridiagonal::SolutionState{
						.u = span<double>(uf.data() + k * nz, nz),
						.d = span<double>(ud.data() + k * nz, nz),
					}
				);

				matrix.solve(
					Tridiagonal::SolutionState{
						.u = span<double>(vf.data() + k * nz, nz),
						.d = span<double>(vd.data() + k * nz, nz),
					}
				);
			}
		}
	}
}

void Solver::calcZ(
	int nx, int ny, int nz, 
	double dt, double dx, double dy, 
	vector<double>& h, 
	vector<double>& ua, 
	vector<double>& va, 
	vector<double>& z
) {
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			size_t k = j + i * ny;

			if (h[k] > 0) {
				if (i > 0 && i < nx - 1) {
					int kr = j + (i + 1) * ny;
					int kl = j + (i - 1) * ny;

					z[k] -= dt * (h[kr] * ua[kr] - h[kl] * ua[kl]) / 2 / dx;
				}

				if (j > 0 && j < ny - 1) {
					int ku = j + 1 + i * ny;
					int kd = j - 1 + i * ny;

					z[k] -= dt * (h[ku] * va[ku] - h[kd] * va[kd]) / 2 / dy;
				}
			}
		}
	}
}

double Solver::maxAbsResidual(
	int nx, int ny, int nz, 
	vector<double>& h,
	vector<double>& al, 
	vector<double>& ac, 
	vector<double>& ar, 
	vector<double>& uf, 
	vector<double>& vf, 
	vector<double>& ud, 
	vector<double>& vd
) {
	double maxRes = 0;

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			int k = j + i * ny;

			if (h[k] > 0) {
				Tridiagonal::Matrix matrix(
					Tridiagonal::MatrixParams{
						.l = span<double>(al.data() + k * nz, nz),
						.c = span<double>(ac.data() + k * nz, nz),
						.r = span<double>(ar.data() + k * nz, nz)
					}
				);

				Tridiagonal::ResidualState res{
					.u = span<double>(uf.data() + k * nz, nz),
					.d = span<double>(ud.data() + k * nz, nz),
				};

				matrix.residual(res);

				for (int k = 0; k < nz; k++) {
					maxRes = max(maxRes, abs(res.d[k]));
				}

				res = Tridiagonal::ResidualState{
					.u = span<double>(vf.data() + k * nz, nz),
					.d = span<double>(vd.data() + k * nz, nz),
				};

				matrix.residual(res);

				for (int k = 0; k < nz; k++) {
					maxRes = max(maxRes, abs(res.d[k]));
				}
			}			
		}
	}

	return maxRes;
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
			Column {.name = "va", .precision = 3 },
			Column {.name = "w10", .precision = 3 },
			Column {.name = "q", .precision = 3 },
			Column {.name = "z", .precision = 7 }
		},
		vector<Writers::Surface::DataKind>{
			SurfaceVector {
				.vx = current.ua,
				.vy = current.va
			},
			SurfaceVector {
				.vx = currentWindData.speed.u10,
				.vy = currentWindData.speed.v10,
			},
			SurfaceVector {
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
			Column {.name = "v", .precision = 3 }
		},
		vector<Writers::Volume::DataKind>{
			VolumeVector {.vx = u, .vy = v, .vz = w }
		},
		dirs.volume
	);

	volumeNodeWriter.write(
		TimeState{ .m = time.m, .t = time.t },
		VolumeState{ .nx = volume.nx, .ny = volume.ny, .nz = volume.nz + 1, .dx = volume.dx, .dy = volume.dy, .dz = volume.dz, .h = volume.h },
		vector<Column>{
			Column {.name = "nu", .precision = 3 }
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
	vector<double> u1a(nx * ny);

	vector<double> va(nx * ny);
	vector<double> v1a(nx * ny);

	vector<double> z(nx * ny);

	vector<double> uf(nx * ny * nz);
	vector<double> ud(nx * ny * nz);
	vector<double> udt(nx * ny * nz);

	vector<double> vf(nx * ny * nz);
	vector<double> vd(nx * ny * nz);
	vector<double> vdt(nx * ny * nz);

	vector<double> u(nx * ny * nz);
	vector<double> v(nx * ny * nz);
	vector<double> w(nx * ny * nz);

	vector<double> al(nx * ny * nz);
	vector<double> ac(nx * ny * nz);
	vector<double> ar(nx * ny * nz);

	setInitialCondition(
		nx, ny, nz,
		uf, vf,
		ua, va,
		z
	);

	auto windData = qGenerator->generate(nx, ny);

	size_t currentWindIndex = 0;

	auto nu = nuGenerator->generate(
		VolumeState{
			.nx = nx, .ny = ny, .nz = nz + 1,
			.dx = dx, .dy = dy, 
			.dz = dz, .h = h
		},
		WindCurrent{
			.speed = WindSpeed {
				.u10 = windData[currentWindIndex].speed.u10,
				.v10 = windData[currentWindIndex].speed.v10
			},
			.stress = WindStress {
				.qx = windData[currentWindIndex].stress.qx,
				.qy = windData[currentWindIndex].stress.qy
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

	vector<double> zp(nx * ny);

	vector<double> up(nx * ny * nz);
	vector<double> vp(nx * ny * nz);

	Utils::Data::updateData(nx, ny, zp, z);
	Utils::Data::updateData(nx, ny, nz, up, uf);
	Utils::Data::updateData(nx, ny, nz, vp, vf);

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
			Column { .name = "h", .precision = 3 }
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

	long long calcTime = 0;

	auto start = high_resolution_clock::now();

	while (t <= endTime) {
		calcUA(
			nx, ny, nz,
			phys, kb, g,
			dx, dt, h,
			windData[currentWindIndex].stress.qx,
			ua, va, z, 
			uf, u1a
		);

		calcVA(
			nx, ny, nz,
			phys, kb, g,
			dy, dt, h,
			windData[currentWindIndex].stress.qy,
			ua, va, z,
			vf, v1a
		);

		calcRHS(
			nx, ny, nz, 
			phys, kb,
			dx, dy, dt,
			dz, h,
			nu, 
			windData[currentWindIndex].stress.qx, windData[currentWindIndex].stress.qy,
			ua, u1a, va, v1a, 
			uf, vf, ud, vd
		);

		udt = ud;
		vdt = vd;

		createTridiagonalMatrix(
			kb, dt, 
			nx, ny, nz,
			dz, h, nu, 
			al, ac, ar
		);

		calcUVF(
			nx, ny, nz, 
			h, 
			al, ac, ar,
			uf, vf,
			ud, vd
		);

		swap(ua, u1a);
		swap(va, v1a);

		calcZ(
			nx, ny, nz, 
			dt, dx, dy,
			h, 
			ua, va, z
		);

		t += dt;
		n += 1;

		if (currentWindIndex < windData.size() - 1 && t >= windData[currentWindIndex + 1].time) {
			currentWindIndex += 1;			

			nu = nuGenerator->generate(
				VolumeState{
					.nx = nx, .ny = ny, .nz = nz + 1,
					.dx = dx, .dy = dy,
					.dz = dz, .h = h
				},
				WindCurrent{
					.speed = WindSpeed {
						.u10 = windData[currentWindIndex].speed.u10,
						.v10 = windData[currentWindIndex].speed.v10
					},
					.stress = WindStress {
						.qx = windData[currentWindIndex].stress.qx,
						.qy = windData[currentWindIndex].stress.qy
					},
					.current = Current {
						.ua = ua,
						.va = va
					}
				}
			);			

			nuMax = Utils::Viscosity::maxNU(nu, nx, ny, nz + 1);

			dt = timeStep->calculate(Time::State{
				.t = t,
				.dt = dtp,
				.dtMax = dtMax,
				.mult = false
			});

			cout << format("Update wind data dt = {}", dt) << endl;
		}

		double dt1 = min(dtMax, hMin * dzMin * hMin * dzMin / 2 / nuMax);

		if (dt1 < dtp) {
			dtp = dt1;

			dt = timeStep->calculate(Time::State{
				.t = t,
				.dt = dt1,
				.dtMax = dtMax,
				.mult = false
			});
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

			auto umd = Utils::Data::maxAbsDifference(nx, ny, nz, up, uf);
			auto vmd = Utils::Data::maxAbsDifference(nx, ny, nz, vp, vf);

			auto maxRes = maxAbsResidual(nx, ny, nz, h, al, ac, ar, uf, vf, udt, vdt);

			cout << format(
				"Write data in file t={:.3f}, nu={:.5f}, maxAbsRes={}, convergence of u={:.5f}, v={:.5f} with dt={:.5}, calc time={}",
				t, nuMax, maxRes, umd, vmd, dt, calcTime / 1000
			) << endl;		

			statistics.push_back({ (double)n, tn, (double)(calcTime / 1000), maxRes, umd, vmd });

			Utils::Data::updateData(nx, ny, zp, z);
			Utils::Data::updateData(nx, ny, nz, up, uf);
			Utils::Data::updateData(nx, ny, nz, vp, vf);

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
			Writers::Statistics::Column{ .name = "vmd", .precision = 5 }
		},
		statistics,
		dirs.root
	);

	cout << "Total number of iterations: " << n;
}
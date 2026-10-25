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
	double dx, double dy, bool symmetrise,
	bool calcNU, bool interpolateWind,
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

	this->symmetrise = symmetrise;

	this->calcNU = calcNU;
	this->interpolateWind = interpolateWind;

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
	size_t nx, size_t ny, size_t nz,
	vector<double>& uf,
	vector<double>& vf,
	vector<double>& ua,
	vector<double>& va,
	vector<double>& z
) {
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			size_t p = j + i * ny;

			for (size_t k = 0; k < nz; k++) {
				size_t id = k + p * nz;

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
	size_t nx, size_t ny, size_t nz,
	PhysicalParams phys, 
	double kb, double g, 
	double dx, double dt, 
	vector<double>& h, 
	vector<double>& dz,
	vector<double>& qx, 
	vector<double>& ua, 
	vector<double>& va, 
	vector<double>& z, 
	vector<double>& uf,
	vector<double>& nu,
	vector<double>& u1a
) {
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			size_t k = j + i * ny;
			size_t idx = nz - 1 + nz * k;

			if (h[k] > 0) {
				size_t nb = nz + k * (nz + 1);

				double ab = kb * h[k] * dz[nz - 1] / 2 / nu[nb];
				double ub = (uf[idx] + ua[k]) / (1 + ab);

				u1a[k] = ua[k] + phys.f * va[k] * dt - kb * dt * ub / h[k] + qx[k] * dt / phys.rho / h[k];

				size_t kr = k;
				size_t kl = k;

				if (i == 0) {
					kr = j + ny;
					kl = j + (nx - 1) * ny;					
				} else if (i == nx - 1) {
					kr = j;
					kl = j + (nx - 2) * ny;
				} else {
					kr = j + (i + 1) * ny;
					kl = j + (i - 1) * ny;
				}

				double gx = (phys.f * va[k] + qx[k] / phys.rho / h[k]) / g;

				double zr = (h[kr] > 0) ? z[kr] : z[k] + dx * gx;
				double zl = (h[kl] > 0) ? z[kl] : z[k] - dx * gx;

				u1a[k] -= g * dt * (zr - zl) / 2 / dx;
			}
		}
	}
}

void Solver::calcVA(
	size_t nx, size_t ny, size_t nz,
	PhysicalParams phys, 
	double kb, double g, 
	double dy, double dt, 
	vector<double>& h, 
	vector<double>& dz,
	vector<double>& qy, 
	vector<double>& ua, 
	vector<double>& va, 
	vector<double>& z, 
	vector<double>& vf, 
	vector<double>& nu,
	vector<double>& v1a
) {
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			size_t k = j + i * ny;
			size_t idx = nz - 1 + nz * k;

			if (h[k] > 0) {
				size_t nb = nz + k * (nz + 1);

				double ab = kb * h[k] * dz[nz - 1] / 2 / nu[nb];
				double vb = (vf[idx] + va[k]) / (1 + ab);

				v1a[k] = va[k] - phys.f * ua[k] * dt - kb * dt * vb / h[k] + qy[k] * dt / phys.rho / h[k];

				size_t ku = k;
				size_t kd = k;

				if (j == 0) {
					ku = 1 + i * ny;
					kd = ny - 1 + i * ny;
				}
				else if (j == ny - 1) {
					ku = i * ny;
					kd = ny - 2 + i * ny;
				} else {
					ku = j + 1 + i * ny;
					kd = j - 1 + i * ny;					
				}

				double gy = (-phys.f * ua[k] + qy[k] / phys.rho / h[k]) / g;

				double zu = (h[ku] > 0) ? z[ku] : z[k] + dy * gy;
				double zd = (h[kd] > 0) ? z[kd] : z[k] - dy * gy;

				v1a[k] -= g * dt * (zu - zd) / 2 / dy;
			}
		}
	}
}

void Solver::calcRHS(
	size_t nx, size_t ny, size_t nz,
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

					size_t idu = k + 1 + p * nz;
					size_t idd = k - 1 + p * nz;

					size_t idb = nz - 1 + p * nz;

					size_t nb = nz + p * (nz + 1);

					double ab = kb * h[p] * dz[nz - 1] / 2 / nu[nb];

					double ub = (uf[idb] + ua[p]) / (1 + ab);
					double vb = (vf[idb] + va[p]) / (1 + ab);

					size_t nid = k + p * (nz + 1);
					size_t nidu = k + 1 + p * (nz + 1);

					if (k == 0) {
						double um1 = uf[id] + h[p] * qx[p] * dz[k] / phys.rho / nu[nid];
						double vm1 = vf[id] + h[p] * qy[p] * dz[k] / phys.rho / nu[nid];

						ud[id] = uf[id] + phys.f * vf[id] * dt + dt * (
							nu[nidu] * (uf[idu] - uf[id]) / (dz[k + 1] + dz[k]) -
							nu[nid] * (uf[id] - um1) / 2 / dz[k]
						) / h[p] / h[p] / dz[k] + kb * ub * dt / h[p] - qx[p] * dt / phys.rho / h[p]
						+ dt * qx[p] / phys.rho / h[p] / 2 / dz[k];

						vd[id] = vf[id] - phys.f * uf[id] * dt + dt * (
							nu[nidu] * (vf[idu] - vf[id]) / (dz[k + 1] + dz[k]) -
							nu[nid] * (vf[id] - vm1) / 2 / dz[k]
						) / h[p] / h[p] / dz[k] + kb * vb * dt / h[p] - qy[p] * dt / phys.rho / h[p]
						+ dt * qy[p] / phys.rho / h[p] / 2 / dz[k];
					}
					else if (k == nz - 1) {
						double un = (
							uf[id] - kb * h[p] * dz[k] * uf[id] / 2 / nu[nidu]
							- kb * h[p] * dz[k] * ua[p] / nu[nidu]
						) / (1 + kb * h[p] * dz[k] / 2 / nu[nidu]);

						double vn = (
							vf[id] - kb * h[p] * dz[k] * vf[id] / 2 / nu[nidu]
							- kb * h[p] * dz[k] * va[p] / nu[nidu]
						) / (1 + kb * h[p] * dz[k] / 2 / nu[nidu]);

						ud[id] = uf[id] + phys.f * vf[id] * dt + dt * (
							nu[nidu] * (un - uf[id]) / 2 / dz[k] -
							nu[nid] * (uf[id] - uf[idd]) / (dz[k] + dz[k - 1])
						) / h[p] / h[p] / dz[k] + kb * ub * dt / h[p] - qx[p] * dt / phys.rho / h[p]
						- dt * kb * ua[p] / h[p] / 2 / dz[k] / (1 + kb * h[p] * dz[k] / 2 / nu[nidu]);

						vd[id] = vf[id] - phys.f * uf[id] * dt + dt * (
							nu[nidu] * (vn - vf[id]) / 2 / dz[k] -
							nu[nid] * (vf[id] - vf[idd]) / (dz[k] + dz[k - 1])
						) / h[p] / h[p] / dz[k] + kb * vb * dt / h[p] - qy[p] * dt / phys.rho / h[p]
						- dt * kb * va[p] / h[p] / 2 / dz[k] / (1 + kb * h[p] * dz[k] / 2 / nu[nidu]);
					}
					else {						
						ud[id] = uf[id] + phys.f * vf[id] * dt + dt * (
							nu[nidu] * (uf[idu] - uf[id]) / (dz[k + 1] + dz[k]) -
							nu[nid] * (uf[id] - uf[idd]) / (dz[k] + dz[k - 1])
						) / h[p] / h[p] / dz[k] + kb * ub * dt / h[p] - qx[p] * dt / phys.rho / h[p];

						vd[id] = vf[id] - phys.f * uf[id] * dt + dt * (
							nu[nidu] * (vf[idu] - vf[id]) / (dz[k + 1] + dz[k]) -
							nu[nid] * (vf[id] - vf[idd]) / (dz[k] + dz[k - 1])
						) / h[p] / h[p] / dz[k] + kb * vb * dt / h[p] - qy[p] * dt / phys.rho / h[p];
					}
				}
			}					
		}
	}
}

void Solver::createTridiagonalMatrix(
	double kb, double dt, 
	size_t nx, size_t ny, size_t nz,
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

					size_t nid = k + p * (nz + 1);
					size_t nidu = k + 1 + p * (nz + 1);

					if (k == 0) {
						ac[id] = 1 + dt * nu[nidu] / h[p] / h[p] / dz[k] / (dz[k + 1] + dz[k]);
						ar[id] = -dt * nu[nidu] / h[p] / h[p] / dz[k] / (dz[k + 1] + dz[k]);
					}
					else if (k == nz - 1) {
						al[id] = -dt * nu[nid] / h[p] / h[p] / dz[k] / (dz[k] + dz[k - 1]);
						ac[id] = 1 + dt * (
							kb * h[p] * dz[k] / (1 + kb * h[p] * dz[k] / 2 / nu[nidu]) / 2 / dz[k] +
							nu[nid] / (dz[k] + dz[k - 1])
						) / h[p] / h[p] / dz[k];
					}
					else {
						al[id] = -dt * nu[nid] / h[p] / h[p] / dz[k] / (dz[k] + dz[k - 1]);
						ac[id] = 1 + dt * (
							nu[nidu] / (dz[k + 1] + dz[k]) +
							nu[nid] / (dz[k] + dz[k - 1])
						) / h[p] / h[p] / dz[k];
						ar[id] = -dt * nu[nidu] / h[p] / h[p] / dz[k] / (dz[k + 1] + dz[k]);
					}
				}
			}
		}
	}
}

void Solver::calcUVF(
	size_t nx, size_t ny, size_t nz,
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
	size_t nx, size_t ny, size_t nz,
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
				size_t kr = j + (i + 1) * ny;
				size_t kl = j + (i - 1) * ny;

				size_t ku = j + 1 + i * ny;
				size_t kd = j - 1 + i * ny;

				if (i == 0) {
					if (j == 0) {
						kr = ny;
						kl = (nx - 1) * ny;

						ku = 1;
						kd = ny - 1;
					}
					else if (j == ny - 1) {
						kr = ny - 1 + ny;
						kl = ny - 1 + (nx - 1) * ny;

						ku = 0;
						kd = ny - 2;
					}
					else {
						kr = j + ny;
						kl = j + (nx - 1) * ny;

						ku = j + 1;
						kd = j - 1;
					}
				}
				else if (i == nx - 1) {
					if (j == 0) {
						kr = 0;
						kl = (nx - 2) * ny;

						ku = 1 + (nx - 1) * ny;
						kd = ny - 1 + (nx - 1) * ny;
					}
					else if (j == ny - 1) {
						kr = ny - 1;
						kl = ny - 1 + (nx - 2) * ny;

						ku = (nx - 1) * ny;
						kd = ny - 2 + (nx - 1) * ny;;
					}
					else {
						kr = j;
						kl = j + (nx - 2) * ny;

						ku = j + 1 + (nx - 1) * ny;
						kd = j - 1 + (nx - 1) * ny;
					}
				}
				else if (j == 0) {
					kr = (i + 1) * ny;
					kl = (i - 1) * ny;

					ku = 1 + i * ny;
					kd = ny - 1 + i * ny;
				}
				else if (j == ny - 1) {
					kr = ny - 1 + (i + 1) * ny;
					kl = ny - 1 + (i - 1) * ny;

					ku = i * ny;
					kd = ny - 2 + i * ny;
				}

				double qr = (h[kr] > 0) ? h[kr] * ua[kr] : -h[k] * ua[k];
				double ql = (h[kl] > 0) ? h[kl] * ua[kl] : -h[k] * ua[k];
				double qu = (h[ku] > 0) ? h[ku] * va[ku] : -h[k] * va[k];
				double qd = (h[kd] > 0) ? h[kd] * va[kd] : -h[k] * va[k];

				z[k] -= dt * (qr - ql) / 2 / dx;
				z[k] -= dt * (qu - qd) / 2 / dy;
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

void Solver::interpolateQUV10(
	int nx, int ny, double alpha,
	const vector<double>& u10c,
	const vector<double>& v10c,
	const vector<double>& qxc,
	const vector<double>& qyc,
	const vector<double>& u10n,
	const vector<double>& v10n,
	const vector<double>& qxn,
	const vector<double>& qyn,
	vector<double>& u10,
	vector<double>& v10,
	vector<double>& qx,
	vector<double>& qy
) {
	Utils::Data::interpolate(nx, ny, alpha, u10c, u10n, u10);
	Utils::Data::interpolate(nx, ny, alpha, v10c, v10n, v10);

	Utils::Data::interpolate(nx, ny, alpha, qxc, qxn, qx);
	Utils::Data::interpolate(nx, ny, alpha, qyc, qyn, qy);
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
	size_t nx = symmetrise ? volume.nx / 2 : volume.nx;
	size_t ny = symmetrise ? volume.ny / 2 : volume.ny;

	surfaceWriter.write(
		TimeState{ .m = time.m, .t = time.t },
		SurfaceState{ .nx = nx, .ny = ny, .bathymetry = volume.bathymetry },
		vector<Column>{
			Column {.name = "va", .precision = 9 },
			Column {.name = "w10", .precision = 9 },
			Column {.name = "q", .precision = 9 },
			Column {.name = "z", .precision = 18 }
		},
		vector<Writers::Surface::DataKind>{
			SurfaceVector {
				.vx = currentWindData.current.ua,
				.vy = currentWindData.current.va
			},
			SurfaceVector {
				.vx = currentWindData.speed.u10,
				.vy = currentWindData.speed.v10,
			},
			SurfaceVector {
				.vx = currentWindData.stress.qx,
				.vy = currentWindData.stress.qy
			},
			Scalar{
				.v = z
			}			
		},
		dirs.surface
	);

	Utils::Data::volumePlusSurface(nx, ny, volume.nz, uf, ua, u);
	Utils::Data::volumePlusSurface(nx, ny, volume.nz, vf, va, v);

	volumeElemWriter.write(
		TimeState{ .m = time.m, .t = time.t },
		VolumeState{ .nx = nx, .ny = ny, .nz = volume.nz, .dz = volume.dz, .bathymetry = volume.bathymetry },
		vector<Column>{
			Column {.name = "v", .precision = 9 }
		},
		vector<Writers::Volume::DataKind>{
			VolumeVector {.vx = u, .vy = v, .vz = w }
		},
		dirs.volume
	);

	volumeNodeWriter.write(
		TimeState{ .m = time.m, .t = time.t },
		VolumeState{ .nx = nx, .ny = ny, .nz = volume.nz + 1, .dz = volume.dz, .bathymetry = volume.bathymetry },
		vector<Column>{
			Column {.name = "nu", .precision = 9 }
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

	auto geo = areaGenerator->generate();

	size_t nx = static_cast<int>(
		ceil(geo.l / dx)
	) + 1;

	size_t ny = static_cast<int>(
		ceil(geo.w / dy)
	) + 1;

	auto dz = dzGenerator->generate();

	double dzMin = *min_element(dz.begin(), dz.end());

	size_t nz = dz.size();

	auto bathymetry = hGenerator->generate(nx, ny, dx, dy);

	auto [hMin, hMax] = Utils::Bathymetry::minMaxH(bathymetry.h, nx, ny);

	auto windData = qGenerator->generate(nx, ny, dx, dy);

	if (symmetrise) {
		bathymetry.h = Utils::Data::symmetrise(nx, ny, bathymetry.h);

		for (auto& wind : windData) {
			wind.speed.u10 = Utils::Data::symmetrise(nx, ny, wind.speed.u10);
			wind.speed.v10 = Utils::Data::symmetrise(nx, ny, wind.speed.v10);

			wind.stress.qx = Utils::Data::symmetrise(nx, ny, wind.stress.qx);
			wind.stress.qy = Utils::Data::symmetrise(nx, ny, wind.stress.qy);
		}

		nx *= 2;
		ny *= 2;
	}	

	const double dtCFL = min(dx, dy) / sqrt(2 * g * hMax) / 1.5;
	const double dtMax = outputTimeStep / std::ceil(outputTimeStep / dtCFL);

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

	size_t currentWindIndex = 0;

	WindSpeed windSpeed{
		.u10 = windData[currentWindIndex].speed.u10,
		.v10 = windData[currentWindIndex].speed.v10
	};

	WindStress windStress{
		.qx = windData[currentWindIndex].stress.qx,
		.qy = windData[currentWindIndex].stress.qy
	};

	Current current{
		.ua = ua,
		.va = va
	};

	auto nu = nuGenerator->generate(
		VolumeState{
			.nx = nx, .ny = ny, .nz = nz + 1,
			.dz = dz,
			.bathymetry = bathymetry
		},
		WindCurrent{
			.speed = windSpeed,
			.stress = windStress,
			.current = current
		}
	);

	double nuMax = Utils::Viscosity::maxNU(nu, nx, ny, nz + 1);

	double dt = min(dtMax, hMin * dzMin * hMin * dzMin / 2 / nuMax);

	bool updateDT = false;

	vector<double> zp(nx * ny);

	vector<double> up(nx * ny * nz);
	vector<double> vp(nx * ny * nz);

	Utils::Data::updateData(nx, ny, zp, z);
	Utils::Data::updateData(nx, ny, nz, up, uf);
	Utils::Data::updateData(nx, ny, nz, vp, vf);

	double t = 0;
	double tn = outputTimeStep;

	double t0 = windData[currentWindIndex].time;
	double t1 = windData[currentWindIndex].time;

	if (interpolateWind && currentWindIndex + 1 < windData.size()) {
		t1 = windData[currentWindIndex + 1].time;
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
		SurfaceState{ .nx = symmetrise ? nx / 2 : nx, .ny = symmetrise ? ny / 2 : ny, .bathymetry = bathymetry }, {
			Column { .name = "h", .precision = 3 }
		}, {
			Scalar{ .v = bathymetry.h }
		},
		dirs.root
	);

	VolumeState volume = { .nx = nx, .ny = ny, .nz = nz, .dz = dz, .bathymetry = bathymetry };

	writeData(
		TimeState{ .m = m, .t = t },
		volume,
		WindCurrent{
			.speed = windSpeed,
			.stress = windStress,
			.current = current
		},
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

	m += 1;

	while (t < endTime - 1e-9) {
		double dtStep = std::min(dt, tn - t);

		calcUA(
			nx, ny, nz,
			phys, kb, g,
			dx, dtStep, bathymetry.h, dz,
			windStress.qx,
			ua, va, z, 
			uf, nu, u1a
		);

		swap(ua, u1a);

		calcVA(
			nx, ny, nz,
			phys, kb, g,
			dy, dtStep, bathymetry.h, dz,
			windStress.qy,
			ua, va, z,
			vf, nu, v1a
		);

		swap(va, v1a);

		Current currentData{ .ua = ua, .va = va };

		calcRHS(
			nx, ny, nz, 
			phys, kb,
			dx, dy, dtStep,
			dz, bathymetry.h,
			nu, 
			windStress.qx, windStress.qy,
			ua, u1a, va, v1a, 
			uf, vf, ud, vd
		);

		udt = ud;
		vdt = vd;

		createTridiagonalMatrix(
			kb, dtStep, 
			nx, ny, nz,
			dz, bathymetry.h, nu, 
			al, ac, ar
		);

		calcUVF(
			nx, ny, nz, 
			bathymetry.h, 
			al, ac, ar,
			uf, vf,
			ud, vd
		);

		calcZ(
			nx, ny, nz, 
			dtStep, dx, dy,
			bathymetry.h, 
			ua, va, z
		);

		t += dtStep;

		if (std::abs(t - tn) < 1e-9) t = tn;

		n += 1;

		if (currentWindIndex + 1 < windData.size() && t >= windData[currentWindIndex + 1].time) {
			currentWindIndex += 1;	

			windSpeed.u10 = windData[currentWindIndex].speed.u10;
			windSpeed.v10 = windData[currentWindIndex].speed.v10;

			windStress.qx = windData[currentWindIndex].stress.qx;
			windStress.qy = windData[currentWindIndex].stress.qy;

			if (interpolateWind && currentWindIndex + 1 < windData.size()) {
				t0 = windData[currentWindIndex].time;
				t1 = windData[currentWindIndex + 1].time;
			}

			updateDT = !interpolateWind;
		}

		if (interpolateWind) {
			double alpha = (t1 == t0) ? 0 : std::clamp((t - t0) / (t1 - t0), 0.0, 1.0);

			size_t wi = currentWindIndex;
			size_t win = currentWindIndex;

			if (windData.size() > 1) {
				if (currentWindIndex + 1 < windData.size()) {
					win = currentWindIndex + 1;
				}
				else {
					wi = currentWindIndex - 1;
					win = currentWindIndex;
				}
			}

			interpolateQUV10(
				nx, ny, alpha, 
				windData[wi].speed.u10, windData[wi].speed.v10,
				windData[wi].stress.qx, windData[wi].stress.qy,
				windData[win].speed.u10, windData[win].speed.v10,
				windData[win].stress.qx, windData[win].stress.qy,
				windSpeed.u10, windSpeed.v10, windStress.qx, windStress.qy
			);
		}

		if (calcNU) {
			nu = nuGenerator->generate(
				VolumeState{
					.nx = nx, .ny = ny, .nz = nz + 1,
					.dz = dz,
					.bathymetry = bathymetry
				},
				WindCurrent{
					.speed = windSpeed,
					.stress = windStress,
					.current = currentData
				}
			);

			nuMax = Utils::Viscosity::maxNU(nu, nx, ny, nz + 1);
		}

		double dt1 = min(dtMax, hMin * dzMin * hMin * dzMin / 2 / nuMax);

		if (updateDT) {
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

			writeData(
				TimeState{ .m = m, .t = t },
				volume,
				WindCurrent{
					.speed = windSpeed,
					.stress = windStress,
					.current = currentData
				},
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
			auto zm = Utils::Data::maxAbsDifference(nx, ny, zp, z);

			auto maxRes = maxAbsResidual(nx, ny, nz, bathymetry.h, al, ac, ar, uf, vf, udt, vdt);

			std::cout << format(
				"Write data in file t={:.3f}, nu={:.5f}, maxAbsRes={}, convergence of u={:.5f}, v={:.5f}, zm={:.19f} with dt={:.5}, calc time={}",
				t, nuMax, maxRes, umd, vmd, zm, dt, calcTime / 1000
			) << endl;		

			statistics.push_back({ (double)n, tn, (double)(calcTime / 1000), maxRes, umd, vmd, zm });

			Utils::Data::updateData(nx, ny, zp, z);
			Utils::Data::updateData(nx, ny, nz, up, uf);
			Utils::Data::updateData(nx, ny, nz, vp, vf);

			m += 1;
			tn += outputTimeStep;

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
			Writers::Statistics::Column{ .name = "zm", .precision = 19 }
		},
		statistics,
		dirs.root
	);

	std::cout << "Total number of iterations: " << n;
}
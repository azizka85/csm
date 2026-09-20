#include "wind_speed_generator.h"

using namespace Generators::Viscosity;

WindSpeedGenerator::WindSpeedGenerator(
	shared_ptr<Turbulence::IFrictionalInfluence> frictionalInfluence,
	shared_ptr<Turbulence::ITidalSurfaceViscosity> tidalSurfaceViscosity,
	shared_ptr<Turbulence::IViscosity> viscosity,
	string dirAttr
) {
	this->frictionalInfluence = frictionalInfluence;
	this->tidalSurfaceViscosity = tidalSurfaceViscosity;
	this->viscosity = viscosity;

	this->dirAttr = move(dirAttr);
}

path WindSpeedGenerator::addDirectory(path outDir) {
	auto attr = dirAttr.empty() ? "wnu" : dirAttr;

	return outDir / viscosity->dirName(attr);
}

vector<double> WindSpeedGenerator::generate(Data::VolumeState volume, Data::WindCurrent windCurrent) {
	auto nx = volume.nx;
	auto ny = volume.ny;
	auto nz = volume.nz;

	auto& dz = volume.dz;
	auto& h = volume.h;

	auto& u10Data = windCurrent.speed.u10;
	auto& v10Data = windCurrent.speed.v10;

	auto& qxData = windCurrent.stress.qx;
	auto& qyData = windCurrent.stress.qy;

	auto& uaData = windCurrent.current.ua;
	auto& vaData = windCurrent.current.va;

	vector<double> nu(nx * ny * nz, 0);

	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			auto p = j + i * ny;

			auto u10 = u10Data[p];
			auto v10 = v10Data[p];

			auto qx = qxData[p];
			auto qy = qyData[p];

			auto ua = uaData[p];
			auto va = vaData[p];

			auto windCurrent = Turbulence::WindCurrentParams {
				.wind = Wind::SpeedVector { 
					.u10 = u10, 
					.v10 = v10 
				},
				.stress = Wind::StressVector { 
					.qx = qx, 
					.qy = qy 
				},
				.current = Data::DepthAveragedVelocity { 
					.ua = ua, 
					.va = va 
				}
			};

			auto fricInf = frictionalInfluence->calculate(windCurrent);

			auto friction = Turbulence::FrictionParams{
				.fricInf = fricInf,
				.windCurrent = windCurrent
			};

			auto visc = tidalSurfaceViscosity->calculate(friction);

			auto viscState = Turbulence::ViscosityParams { 
				.viscosity = visc,
				.friction = friction
			};

			double z = 0;

			for (size_t k = 0; k < nz; k++) {
				auto id = k + p * nz;

				nu[id] = viscosity->calculate(z, viscState);

				if (k < nz - 1) {
					z += h[p] * dz[k];
				}
			}
		}
	}

	return nu;
}

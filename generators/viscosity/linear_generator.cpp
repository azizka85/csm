#include <utils/checks/phys.h>
#include <utils/checks/bathymetry.h>

#include "linear_generator.h"

using namespace Utils::Checks;

using namespace Generators::Viscosity;

LinearGenerator::LinearGenerator(double ht, Turbulence::ViscosityState viscosity) {
	Bathymetry::checkH(ht);

	Phys::checkNU(viscosity.nus);
	Phys::checkNU(viscosity.nut);

	this->ht = ht;

	this->viscosity = move(viscosity);
}

path LinearGenerator::addDirectory(path outDir) {
	return outDir / path(
		format("lnu, ht={}, nus={}, nut={}", ht, viscosity.nus, viscosity.nut)
	);
}

vector<double> LinearGenerator::generate(Data::VolumeState volume, Data::WindCurrent windCurrent) {
	auto nx = volume.nx;
	auto ny = volume.ny;
	auto nz = volume.nz;

	auto& dz = volume.dz;
	auto& h = volume.h;

	auto nus = viscosity.nus;
	auto nut = viscosity.nut;

	vector<double> nu(nx * ny * nz, 0);

	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			auto p = j + i * ny;

			double z = 0;

			for (size_t k = 0; k < nz; k++) {
				auto id = k + p * nz;

				if (z < ht) {
					nu[id] = nus - (nus - nut) * z / ht;
				}
				else {
					nu[id] = nut;
				}

				if (k < nz - 1) {
					z += h[p] * dz[k];
				}
			}
		}
	}

	return nu;
}

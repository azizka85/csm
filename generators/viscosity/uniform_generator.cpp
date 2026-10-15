#include <utils/checks/phys.h>

#include "uniform_generator.h"

using namespace Utils::Checks;

using namespace Generators::Viscosity;

UniformGenerator::UniformGenerator(double nu) {
	Phys::checkNU(nu);

	this->nu = nu;
}

path UniformGenerator::addDirectory(path outDir) {
	return outDir / path(
		format("unu, nu={}", nu)
	);
}

vector<double> UniformGenerator::generate(const Data::VolumeState& volume, const Data::WindCurrent& windCurrent) {
	auto nx = volume.nx;
	auto ny = volume.ny;
	auto nz = volume.nz;
	
	return vector<double>(nx * ny * nz, nu);
}
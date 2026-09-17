#include <utils/checks/bathymetry.h>

#include "uniform_generator.h"

using namespace Utils::Checks;

using namespace Generators;

DZ::UniformGenerator::UniformGenerator(double dz) {
	Bathymetry::checkDZ(dz);

	this->dz = dz;
}

path DZ::UniformGenerator::addDirectory(path outDir) {
	return outDir / path(
		format("udz, dz={}", dz)
	);
}

vector<double> DZ::UniformGenerator::generate() {
	size_t nz = ceil(1 / dz);

	return vector<double>(nz, dz);
}
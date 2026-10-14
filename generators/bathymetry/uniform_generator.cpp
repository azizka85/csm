#include <utils/checks/bathymetry.h>

#include "uniform_generator.h"

using namespace Utils;

using namespace Generators;

Bathymetry::UniformGenerator::UniformGenerator(double h) {
	Checks::Bathymetry::checkH(h);

	this->h = h;
}

path Bathymetry::UniformGenerator::addDirectory(path outDir) {
	return outDir / path(
		format("uh, h={}", h)
	);
}

vector<double> Bathymetry::UniformGenerator::generate(size_t nx, size_t ny) {
	return vector<double>(nx * ny, h);
}
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

Data::Bathymetry Bathymetry::UniformGenerator::generate(
	size_t nx, size_t ny,
	double dx, double dy
) {
	Data::Bathymetry bathymetry{
		.h = vector<double>(nx * ny, h),
		.x = vector<double>(nx * ny),
		.y = vector<double>(nx * ny)
	};

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			double x = dx * i;
			double y = dy * j;

			size_t p = j + i * ny;

			bathymetry.x[p] = x;
			bathymetry.y[p] = y;
		}
	}

	return bathymetry;
}
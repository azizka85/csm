#define _USE_MATH_DEFINES

#include <cmath>

#include <utils/checks/bathymetry.h>

#include "cosine_generator.h"

using namespace Generators::Bathymetry;

CosineGenerator::CosineGenerator(double hm) {
	Utils::Checks::Bathymetry::checkH(hm);

	this->hm = hm;
}

path CosineGenerator::addDirectory(path outDir) {
	return outDir / path(
		format("ch, h={}", hm)
	);
}

Data::Bathymetry Generators::Bathymetry::CosineGenerator::generate(
    size_t nx, size_t ny,
    double dx, double dy
) {
    Data::Bathymetry bathymetry{
        .h = vector<double>(nx * ny),
        .x = vector<double>(nx * ny),
        .y = vector<double>(nx * ny)
    };

    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            size_t id = j + i * ny;

            double x = dx * i;
            double y = dy * j;

            bathymetry.h[id] = hm * (1 - 0.98 * cos(2 * M_PI * (i / (nx - 1.) + j / (ny - 1.)))) / 2;

            bathymetry.x[id] = x;
            bathymetry.y[id] = y;
        }
    }

    return bathymetry;
}

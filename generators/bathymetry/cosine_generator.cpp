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

vector<double> Generators::Bathymetry::CosineGenerator::generate(double nx, double ny) {
    vector<double> h(nx * ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int id = j + i * ny;

            h[id] = hm * (1 - 0.8 * cos(2 * M_PI * (i / (nx - 1.) + j / (ny - 1.)))) / 2;
        }
    }

    return h;
}

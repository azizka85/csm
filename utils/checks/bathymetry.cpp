#include <format>

#include <stdexcept>

#include "bathymetry.h"

using namespace std;

using namespace Utils::Checks;

void Bathymetry::checkH(double h) {
    if (h <= 0) {
        throw runtime_error(
            format("h should be > 0, but it is {}", h)
        );
    }
}

void Bathymetry::checkMinDepth(double minDepth) {
    if (minDepth < 0) {
        throw runtime_error(
            format("minDepth should be >= 0, but it is {}", minDepth)
        );
    }
}

void Bathymetry::checkDZ(double dz) {
    if (dz <= 0) {
        throw runtime_error(
            format("dz should be > 0, but it is {}", dz)
        );
    }
}

void Bathymetry::checkDZMinMax(double dzMin, double dzMax) {
    checkDZ(dzMin);
    checkDZ(dzMax);

    if (dzMin >= dzMax) {
        throw runtime_error(
            format("dzMin should be < dzMax, but they are {} and {}", dzMin, dzMax)
        );
    }
}

void Bathymetry::checkZ(double z) {
    if (z < 0 || z > 1) {
        throw runtime_error(
            format("z should be in the range [0, 1], but it is {}", z)
        );
    }
}

void Bathymetry::checkZM(double zm) {
    if (zm < 0.5 || zm > 1) {
        throw runtime_error(
            format("zm should be in the range [0.5, 1], but it is {}", zm)
        );
    }
}
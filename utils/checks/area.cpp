#include <format>

#include <stdexcept>

#include "area.h"

using namespace std;

using namespace Utils::Checks;

void Area::checkL(double l) {
    if (l <= 0) {
        throw runtime_error(
            format("l should be > 0, but it is {}", l)
        );
    }
}

void Area::checkW(double w) {
    if (w <= 0) {
        throw runtime_error(
            format("w should be > 0, but it is {}", w)
        );
    }
}

void Area::checkLonMinMax(double lonMin, double lonMax) {
    if (lonMin >= lonMax) {
        throw runtime_error(
            format("lonMin should be < lonMax, but they are {} and {}", lonMin, lonMax)
        );
    }
}

void Area::checkLatMinMax(double latMin, double latMax) {
    if (latMin >= latMax) {
        throw runtime_error(
            format("latMin should be < latMax, but they are {} and {}", latMin, latMax)
        );
    }

}

void Area::checkNumNodes(size_t nx, size_t ny) {
    if (nx == 0 || ny == 0) {
        throw runtime_error(
            format("nx and ny should be > 0, but they are {} and {}", nx, ny)
        );
    }
}

void Area::checkNumElemNodes(size_t nx, size_t ny) {
    if (nx * ny <= 1) {
        throw runtime_error(
            format("number of element nodes nx * ny should be > 1, but it is {}", nx * ny)
        );
    }
}

#include <format>

#include <stdexcept>

#include "phys.h"

using namespace std;

using namespace Utils::Checks;

void Phys::checkRho(double rho) {
    if (rho <= 0) {
        throw runtime_error(
            format("rho should be > 0, but it is {}", rho)
        );
    }
}

void Phys::checkCD(double cd) {
    if (cd <= 0) {
        throw runtime_error(
            format("cd should be > 0, but it is {}", cd)
        );
    }
}

void Phys::checkNU(double nu) {
    if (nu <= 0) {
        throw runtime_error(
            format("nu should be > 0, but it is {}", nu)
        );
    }
}

void Phys::checkK0(double k0) {
    if (k0 <= 0) {
        throw runtime_error(
            format("k0 should be > 0, but it is {}", k0)
        );
    }
}

void Phys::checkK(double k) {
    if (k <= 0) {
        throw runtime_error(
            format("k should be > 0, but it is {}", k)
        );
    }
}

void Phys::checkSigma(double sigma) {
    if (sigma <= 0) {
        throw runtime_error(
            format("sigma should be > 0, but it is {}", sigma)
        );
    }
}

void Phys::checkF(double f) {
    if (f <= 0) {
        throw runtime_error(
            format("f should be > 0, but it is {}", f)
        );
    }
}
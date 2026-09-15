#include <format>

#include <stdexcept>

#include "time.h"

using namespace std;

using namespace Utils::Checks;

void Time::checkB(double b) {
    if (b <= 0) {
        throw runtime_error(
            format("b should be > 0, but it is {}", b)
        );
    }
}

void Time::checkTMax(double tMax) {
    if (tMax <= 0) {
        throw runtime_error(
            format("tMax should be > 0, but it is {}", tMax)
        );
    }
}

void Time::checkDTMax(double dtMax) {
    if (dtMax <= 0) {
        throw runtime_error(
            format("dtMax should be > 0, but it is {}", dtMax)
        );
    }
}

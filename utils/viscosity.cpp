#include <algorithm>
#include <cmath>

#include "viscosity.h"

using namespace Utils;

double Viscosity::maxNU(vector<double>& nu) {
    double nuMax = nu[0];

	const int np = static_cast<int>(nu.size());

    #pragma omp parallel for schedule(static) reduction(max:nuMax)
    for (int k = 0; k < np; k++) {
        nuMax = max(nuMax, nu[k]);
    }

    return nuMax;
}

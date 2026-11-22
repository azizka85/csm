#include <algorithm>
#include <cmath>

#include "data.h"

using namespace Utils;

void Data::shiftIndex(vector<size_t>& v, size_t& i, size_t m) {
    while (i < v.size() && v[i] < m) {
        i += 1;
    }
}

double Data::maxAbsDifference(vector<double>& u, vector<double>& u1) {
    double maxDiff = 0.;

    const int np = static_cast<int>(u.size());

    #pragma omp parallel for schedule(static) reduction(max:maxDiff)
    for (int p = 0; p < np; p++) {
        maxDiff = max(maxDiff, abs(u[p] - u1[p]));
    }

    return maxDiff;
}

double Data::maxAbsElem(vector<double>& u) {
    double maxElem = 0;

    const int np = static_cast<int>(u.size());

    #pragma omp parallel for schedule(static) reduction(max:maxElem)
    for (int p = 0; p < np; p++) {
        maxElem = max(maxElem, abs(u[p]));
    }


    return maxElem;
}

void Data::parallelUpdateData(vector<double>& u, vector<double>& u1) {
    const int n = static_cast<int>(u.size());

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        u1[i] = u[i];
    }
}

void Data::volumePlusSurface(
    int nz,
    vector<double>& u, vector<double>& us,
    vector<double>& v
) {
    const int np = static_cast<int>(us.size());

    #pragma omp parallel for schedule(static)
    for (int p = 0; p < np; p++) {
        for (int k = 0; k < nz; k++) {
            const int id = k + p * nz;

            v[id] = u[id] + us[p];
        }
    }

}

vector<double> Data::symmetrise(int nx, int ny, vector<double> &v) {
    vector<double> vt(4 * nx * ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int pu = ny + j + i * ny;
            int pd = ny - 1 - j + i * ny;

            vt[pu] = v[pd];

            int pl = j + i * ny;
            
            int pdr = j + (2*nx - 1 - i) * ny;
            int pur = ny + j + (2 * nx - 1 - i) * ny;

            vt[pdr] = vt[pl];
            vt[pur] = vt[pu];
        }
    }

    return vt;
}

void Data::interpolate(double alpha, const vector<double>& vc, const vector<double>& vn, vector<double>& v) {
    const int np = static_cast<int>(vc.size());
    
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < np; k++) {
        v[k] = vc[k] + alpha * (vn[k] - vc[k]);
    }

}

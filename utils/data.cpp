#include "data.h"

using namespace Utils;

double Data::maxAbsDifference(int nx, int ny, vector<double>& u, vector<double>& u1) {
    double maxDiff = 0.;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            maxDiff = max(
                maxDiff,
                abs(u[p] - u1[p])
            );
        }
    }

    return maxDiff;
}

double Data::maxAbsDifference(int nx, int ny, int nz, vector<double>& u, vector<double>& u1) {
    double maxDiff = 0.;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            for (int k = 0; k < nz; k++) {
                int id = k + p * nz;

                maxDiff = max(
                    maxDiff,
                    abs(u[id] - u1[id])
                );
            }
        }
    }

    return maxDiff;
}

void Data::updateData(int nx, int ny, vector<double>& u, vector<double>& u1) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            u[p] = u1[p];
        }
    }
}

void Data::updateData(int nx, int ny, int nz, vector<double>& u, vector<double>& u1) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            for (int k = 0; k < nz; k++) {
                int id = k + p * nz;

                u[id] = u1[id];
            }
        }
    }
}

void Data::volumePlusSurface(
    int nx, int ny, int nz,
    vector<double>& u, vector<double>& us,
    vector<double>& v
) {
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int p = j + i * ny;
                int id = k + p * nz;

                v[id] = u[id] + us[p];
            }
        }
    }
}
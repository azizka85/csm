#ifndef UTILS_DATA_H
#define UTILS_DATA_H

#include <vector>

using namespace std;

namespace Utils::Data {
    double maxAbsDifference(
        int nx, int ny,
        vector<double>& u,
        vector<double>& u1
    );

    double maxAbsDifference(
        int nx, int ny, int nz,
        vector<double>& u,
        vector<double>& u1
    );

    void updateData(
        int nx, int ny,
        vector<double>& u,
        vector<double>& u1
    );

    void updateData(
        int nx, int ny, int nz,
        vector<double>& u,
        vector<double>& u1
    );

    void volumePlusSurface(
        int nx, int ny, int nz,
        vector<double>& u, vector<double>& us,
        vector<double>& v
    );
}

#endif
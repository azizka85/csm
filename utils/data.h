#ifndef UTILS_DATA_H
#define UTILS_DATA_H

#include <vector>

using namespace std;

namespace Utils::Data {
    void shiftIndex(vector<size_t>& v, size_t& i, size_t m);

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

    double maxAbsElem(
        int nx, int ny, int nz,
        vector<double>& u
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

    vector<double> symmetrise(int nx, int ny, vector<double> &v);

    void interpolate(int nx, int ny, double alpha, const vector<double>& vc, const vector<double>& vn, vector<double>& v);
}

#endif
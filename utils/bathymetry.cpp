#include "bathymetry.h"

using namespace Utils;

tuple<double, double> Bathymetry::minMaxH(vector<double>& h, int nx, int ny) {
    double minH = 1000000;
    double maxH = h[0];

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int id = j + i * ny;

            if (h[id] > 0) {
                minH = min(minH, h[id]);
                maxH = max(maxH, h[id]);
            }
        }
    }

    return { minH, maxH };
}
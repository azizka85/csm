#ifndef UTILS_BATHYMETRY_H
#define UTILS_BATHYMETRY_H

#include <tuple>
#include <vector>

using namespace std;

namespace Utils::Bathymetry {
	tuple<double, double> minMaxH(vector<double>& h, int nx, int ny);
}

#endif
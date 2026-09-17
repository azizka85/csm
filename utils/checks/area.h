#ifndef UTILS_CHECKS_AREA_H
#define UTILS_CHECKS_AREA_H

namespace Utils::Checks::Area {
	void checkL(double l);
	void checkW(double w);

	void checkLonMinMax(double lonMin, double lonMax);
	void checkLatMinMax(double latMin, double latMax);

	void checkNumNodes(size_t nx, size_t ny);
	void checkNumElemNodes(size_t nx, size_t ny);
}

#endif
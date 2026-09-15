#ifndef UTILS_CHECKS_BATHYMETRY_H
#define UTILS_CHECKS_BATHYMETRY_H

namespace Utils::Checks::Bathymetry {
	void checkH(double h);
	void checkMinDepth(double minDepth);
	
	void checkDZ(double dz);
	void checkDZMinMax(double dzMin, double dzMax);

	void checkZ(double z);
	void checkZM(double zm);
}

#endif
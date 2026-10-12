inline void calc_viscosity(
	const double k, const double sigma, const double num,
	const double u10, const double v10,
	const double ua, const double va,
	double* nus, double* nut
) {
	*nut = k * (ua * ua + va * va) / sigma + num;

	double wind_mag_sq = u10 * u10 + v10 * v10;

	*nus = 0.1825e-3 * pow(wind_mag_sq, 1.25) + num;
}
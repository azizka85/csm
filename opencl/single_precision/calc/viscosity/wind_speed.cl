inline void calc_viscosity(
	const float k, const float sigma, const float num,
	const float u10, const float v10,
	const float ua, const float va,
	float* nus, float* nut
) {
	*nut = k * (ua * ua + va * va) / sigma + num;

	float wind_mag_sq = u10 * u10 + v10 * v10;

	*nus = 0.1825e-3 * pow(wind_mag_sq, 1.25f) + num;
}
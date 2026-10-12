inline void calc_turbulence_frictional_influence(
	const float k0, const  float rho, const float f,
	const float qx, const float qy,
	float* ut, float* ht
) {
	float q = sqrt(qx * qx + qy * qy);

	*ut = sqrt(q / rho);
	*ht = k0 * (*ut) / f;
}

inline void calc_turbulence_frictional_influence(
	const double k0, const  double rho, const double f,
	const double qx, const double qy,
	double* ut, double* ht
) {
	double q = sqrt(qx * qx + qy * qy);

	*ut = sqrt(q / rho);
	*ht = k0 * ut / f;
}

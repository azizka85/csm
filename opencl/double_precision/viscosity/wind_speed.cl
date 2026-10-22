__kernel void viscosity(
	const double k0, const  double rho, const double f,
	const double k, const double sigma, const double num,
	int ny, int nz,
	__global const double* dz,
    __global const double* h,
	__global const double* u10,
    __global const double* v10,
	__global const double* qx,
    __global const double* qy,
	__global const double* ua,
    __global const double* va,
	__global double* nu
) {
	int i = get_global_id(0);
    int j = get_global_id(1);

    int p = j + i*ny;

	if (h[p] > 0) {
		double ut;
		double ht;

		double nus;
		double nut;

		calc_turbulence_frictional_influence(k0, rho, f, qx[p], qy[p], &ut, &ht);
		calc_viscosity(k, sigma, num, u10[p], v10[p], ua[p], va[p], &nus, &nut);

		double z = 0;

		for (int l = 0; l < nz + 1; l++) {
			int id = l + p * (nz + 1);

			nu[id] = (z < ht) ? nus - (nus - nut) * z / ht : nut;

			if (l < nz) {
				z += h[p] * dz[l];
			}
		}
	}
}
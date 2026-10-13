__kernel void viscosity(
	const float k0, const  float rho, const float f,
	const float k, const float sigma, const float num,
	int ny, int nz,
	__global const float* dz,
    __global const float* h,
	__global const float* u10,
    __global const float* v10,
	__global const float* qx,
    __global const float* qy,
	__global const float* ua,
    __global const float* va,
	__global float* nu
) {
	int i = get_global_id(0);
    int j = get_global_id(1);

    int p = j + i*ny;

	if (h[p] > 0) {
		float ut;
		float ht;

		float nus;
		float nut;

		calc_turbulence_frictional_influence(k0, rho, f, qx[p], qy[p], &ut, &ht);
		calc_viscosity(k, sigma, num, u10[p], v10[p], ua[p], va[p], &nus, &nut);

		float z = 0;

		for (int l = 0; l < nz + 1; l++) {
			int id = l + p * (nz + 1);

			nu[id] = (z < ht) ? nus - (nus - nut) * z / ht : nut;

			if (l < nz) {
				z += h[p] * dz[l];
			}
		}
	}
}
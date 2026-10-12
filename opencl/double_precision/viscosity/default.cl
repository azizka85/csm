inline void viscosity(
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

	double ut;
	double ht;

	double nus;
	double nut;

	double z = 0;

	for (int k = 0; k < nz; k++) {
		int id = k + p * nz;

		if (z < ht) {
			nu[id] = nus - (nus - nut) * z / ht;
		} else {
			nu[id] = nut;
		}
	}
}
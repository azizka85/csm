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
) { }
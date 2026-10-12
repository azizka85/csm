inline void viscosity(
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
) { }
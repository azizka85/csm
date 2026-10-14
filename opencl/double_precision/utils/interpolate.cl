__kernel void utils_interpolate_plane(
	int ny,
	const double alpha,	
	__global const double* vc,
	__global const double* vn,
	__global double* v
) {
	int i = get_global_id(0);
    int j = get_global_id(1);

	int k = j + i*ny;

	v[k] = vc[k] + alpha * (vn[k] - vc[k]);
}
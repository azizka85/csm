__kernel void utils_interpolate_plane(
	int ny,
	const float alpha,	
	__global const float* vc,
	__global const float* vn,
	__global float* v
) {
	int i = get_global_id(0);
    int j = get_global_id(1);

	int k = j + i*ny;

	v[k] = vc[k] + alpha * (vn[k] - vc[k]);
}
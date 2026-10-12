__kernel void utils_max_plane(
    int ny, int nz,
     __global const double* u,
      __global double* um
) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    
    int idx = j + i*ny;

    __global const double* us = u + idx*nz;

    double m = us[0];

    for (int k = 1; k < nz; k++) {
        m = fmax(m, us[k]);
    }

    um[idx] = m;
}

__kernel void utils_max_row(
    int ny,
     __global const double* u,
      __global double* um
) {
    int i = get_global_id(0);

    __global const double* us = u + i*ny;

    double m = us[0];

    for (int j = 1; j < ny; j++) {
        m = fmax(m, us[j]);
    }

    um[i] = m;
}

__kernel void utils_max(
    int nx,
     __global const double* u,
      __global double* um
) {
    double m = u[0];

    for (int i = 1; i < nx; i++) {
        m = fmax(m, u[i]);
    }

    um[0] = m;
}

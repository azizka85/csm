inline void slae_direct_tridiagonal_calc(
    int nz,     
    __global const double* l,
    __global const double* c,
    __global const double* r,
    __global double* d,
    __global double* u
) {
    if (nz > 0) {
        u[0] = r[0] / c[0];
        d[0] = d[0] / c[0];

        for (int k = 1; k < nz; k++) {
            double c1 = c[k] - l[k]*u[k-1];
    
            u[k] = r[k] / c1;
            d[k] = (d[k] - l[k]*d[k-1]) / c1;
        }

        u[nz-1] = d[nz-1];

        for (int k = nz-2; k >= 0; k--) {
            u[k] = d[k] - u[k]*u[k+1];
        }
    }
}

inline void slae_direct_tridiagonal_residual(
    int nz,     
    __global const double* l,
    __global const double* c,
    __global const double* r,
    __global const double* u,
    __global double* d
) {
    if (nz > 0) {
        if (nz == 1) {
            d[0] -= c[0] * u[0];
        }
        else {
            d[0] -= c[0] * u[0] + r[0] * u[1];

            for (int k = 1; k < nz - 1; k++) {
                d[k] -= c[k] * u[k] + l[k] * u[k - 1] + r[k] * u[k + 1];
            }

            d[nz - 1] -= c[nz - 1] * u[nz - 1] + l[nz - 1] * u[nz - 2];
        }
    }
}
#include <opencl/single_precision/slae/direct/tridiagonal.cl>

__kernel void wind_induced_currents_davies85_calc_ua(
    float f, float kb, float g, float rho,
    float dx, float dt,
    int nx, int ny, int nz,
    __global const float* h,
    __global const float* qx,
    __global const float* ua,
    __global const float* va,
    __global const float* z,
    __global const float* uf,
    __global float* u1a
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;    
    int idx = nz - 1 + nz*k;

    if (h[k] > 0) {
		u1a[k] = ua[k] + f * va[k] * dt - kb * dt * (uf[idx] + ua[k]) / h[k] + qx[k] * dt / rho / h[k];

		int kr = k;
		int kl = k;

		if (i == 0) {
			kr = j + ny;
			kl = j + (nx - 1) * ny;					
		} else if (i == nx - 1) {
			kr = j;
			kl = j + (nx - 2) * ny;
		} else {
			kr = j + (i + 1) * ny;
			kl = j + (i - 1) * ny;
		}

		u1a[k] -= g * dt * (z[kr] - z[kl]) / 2 / dx;
	}
}

__kernel void wind_induced_currents_davies85_calc_va(
    float f, float kb, float g, float rho,
    float dy, float dt,
    int nx, int ny, int nz,
    __global const float* h,
    __global const float* qy,
    __global const float* ua,
    __global const float* va,
    __global const float* z,
    __global const float* vf,
    __global float* v1a
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;    
    int idx = nz - 1 + nz*k;

    if (h[k] > 0) {
		v1a[k] = va[k] - f * ua[k] * dt - kb * dt * (vf[idx] + va[k]) / h[k] + qy[k] * dt / rho / h[k];

		int ku = k;
		int kd = k;

		if (j == 0) {
			ku = 1 + i * ny;
			kd = ny - 1 + i * ny;
		}
		else if (j == ny - 1) {
			ku = i * ny;
			kd = ny - 2 + i * ny;
		} else {
			ku = j + 1 + i * ny;
			kd = j - 1 + i * ny;					
		}

		v1a[k] -= g * dt * (z[ku] - z[kd]) / 2 / dy;
	}
}

__kernel void wind_induced_currents_davies85_calc_rhs(
    float f, float kb, float rho, 
    float dx, float dy, float dt,
    int ny, int nz,
    __global const float* nu,
    __global const float* dz,
    __global const float* h,
    __global const float* qx,
    __global const float* qy,
    __global const float* ua,
    __global const float* va,
    __global const float* uf,
    __global const float* vf,
    __global float* ud,
    __global float* vd
) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    int k = get_global_id(2);

    int p = j + i*ny;

    if (h[p] > 0) {
        int id = k + p * nz;

		int idu = k + 1 + p * nz;
		int idd = k - 1 + p * nz;

		int idb = nz - 1 + p * nz;

		if (k == 0) {
			float um1 = uf[id] + h[p] * qx[p] * dz[k] / rho / nu[id];
			float vm1 = vf[id] + h[p] * qy[p] * dz[k] / rho / nu[id];

			ud[id] = uf[id] + f * vf[id] * dt + dt * (
				nu[idu] * (uf[idu] - uf[id]) / (dz[k + 1] + dz[k]) -
				nu[id] * (uf[id] - um1) / 2 / dz[k]
			) / h[p] / h[p] / dz[k] + kb * (uf[idb] + ua[p]) * dt / h[p] - qx[p] * dt / rho / h[p]
			+ dt * qx[p] / rho / h[p] / 2 / dz[k];

			vd[id] = vf[id] - f * uf[id] * dt + dt * (
				nu[idu] * (vf[idu] - vf[id]) / (dz[k + 1] + dz[k]) -
				nu[id] * (vf[id] - vm1) / 2 / dz[k]
			) / h[p] / h[p] / dz[k] + kb * (vf[idb] + va[p]) * dt / h[p] - qy[p] * dt / rho / h[p]
			+ dt * qy[p] / rho / h[p] / 2 / dz[k];
		}
		else if (k == nz - 1) {
			float un = (
				uf[id] - kb * h[p] * dz[k] * uf[id] / 2 / nu[idu]
				- kb * h[p] * dz[k] * ua[p] / nu[idu]
			) / (1 + kb * h[p] * dz[k] / 2 / nu[idu]);

			float vn = (
				vf[id] - kb * h[p] * dz[k] * vf[id] / 2 / nu[idu]
				- kb * h[p] * dz[k] * va[p] / nu[idu]
			) / (1 + kb * h[p] * dz[k] / 2 / nu[idu]);

			ud[id] = uf[id] + f * vf[id] * dt + dt * (
				nu[idu] * (un - uf[id]) / 2 / dz[k] -
				nu[id] * (uf[id] - uf[idd]) / (dz[k] + dz[k - 1])
			) / h[p] / h[p] / dz[k] + kb * (uf[idb] + ua[p]) * dt / h[p] - qx[p] * dt / rho / h[p]
			- dt * kb * ua[p] / h[p] / 2 / dz[k] / (1 + kb * h[p] * dz[k] / 2 / nu[idu]);

			vd[id] = vf[id] - f * uf[id] * dt + dt * (
				nu[idu] * (vn - vf[id]) / 2 / dz[k] -
				nu[id] * (vf[id] - vf[idd]) / (dz[k] + dz[k - 1])
			) / h[p] / h[p] / dz[k] + kb * (vf[idb] + va[p]) * dt / h[p] - qy[p] * dt / rho / h[p]
			- dt * kb * va[p] / h[p] / 2 / dz[k] / (1 + kb * h[p] * dz[k] / 2 / nu[idu]);
		}
		else {
			float u = uf[idb] + ua[p];
			float v = vf[idb] + va[p];

			ud[id] = uf[id] + f * vf[id] * dt + dt * (
				nu[idu] * (uf[idu] - uf[id]) / (dz[k + 1] + dz[k]) -
				nu[id] * (uf[id] - uf[idd]) / (dz[k] + dz[k - 1])
			) / h[p] / h[p] / dz[k] + kb * u * dt / h[p] - qx[p] * dt / rho / h[p];

			vd[id] = vf[id] - f * uf[id] * dt + dt * (
				nu[idu] * (vf[idu] - vf[id]) / (dz[k + 1] + dz[k]) -
				nu[id] * (vf[id] - vf[idd]) / (dz[k] + dz[k - 1])
			) / h[p] / h[p] / dz[k] + kb * v * dt / h[p] - qy[p] * dt / rho / h[p];
		}        
    }    
}

__kernel void wind_induced_currents_davies85_create_tridiagonal_matrix(
    float kb, float dt,
    int ny, int nz,
    __global const float* nu,
    __global const float* dz,
    __global const float* h,
    __global float* al,
    __global float* ac,
    __global float* ar
) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    int k = get_global_id(2);

    int p = j + i*ny;

    if (h[p] > 0) {
		int id = k + p * nz;
		int idu = k + 1 + p * nz;

		if (k == 0) {
			ac[id] = 1 + dt * nu[idu] / h[p] / h[p] / dz[k] / (dz[k + 1] + dz[k]);
			ar[id] = -dt * nu[idu] / h[p] / h[p] / dz[k] / (dz[k + 1] + dz[k]);
		}
		else if (k == nz - 1) {
			al[id] = -dt * nu[id] / h[p] / h[p] / dz[k] / (dz[k] + dz[k - 1]);
			ac[id] = 1 + dt * (
				kb * h[p] * dz[k] / (1 + kb * h[p] * dz[k] / 2 / nu[idu]) / 2 / dz[k] +
				nu[id] / (dz[k] + dz[k - 1])
			) / h[p] / h[p] / dz[k];
		}
		else {
			al[id] = -dt * nu[id] / h[p] / h[p] / dz[k] / (dz[k] + dz[k - 1]);
			ac[id] = 1 + dt * (
				nu[idu] / (dz[k + 1] + dz[k]) +
				nu[id] / (dz[k] + dz[k - 1])
			) / h[p] / h[p] / dz[k];
			ar[id] = -dt * nu[idu] / h[p] / h[p] / dz[k] / (dz[k + 1] + dz[k]);
		}
	}
}

__kernel void wind_induced_currents_davies85_calc_uvf(
    int ny, int nz,
	__global const float* h,
    __global const float* al,
    __global const float* ac,
    __global const float* ar,
    __global float* uf,
    __global float* vf,
    __global float* ud,
    __global float* vd
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;

    if (h[k] > 0) {
		int shift = nz * k;

        slae_direct_tridiagonal_calc(
            nz,             
            al + shift, ac + shift, ar + shift, 
            ud + shift, uf + shift
        );

        slae_direct_tridiagonal_calc(
            nz,             
            al + shift, ac + shift, ar + shift, 
            vd + shift, vf + shift
        );        
    }    
}

__kernel void wind_induced_currents_davies85_calc_residual(
	int ny, int nz,
	__global const float* h,
    __global const float* al,
    __global const float* ac,
    __global const float* ar,
    __global const float* uf,
    __global const float* vf,
    __global float* ud,
    __global float* vd
) {
	int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;

    if (h[k] > 0) {
		int shift = nz * k;

        slae_direct_tridiagonal_residual(
            nz,             
            al + shift, ac + shift, ar + shift, 
            uf + shift, ud + shift
        );

        slae_direct_tridiagonal_residual(
            nz,             
            al + shift, ac + shift, ar + shift, 
            vf + shift, vd + shift
        );        
    } 
}

__kernel void wind_induced_currents_davies85_calc_z(
    float dx, float dy, float dt,
    int nx, int ny, int nz,
    __global const float* h,
    __global const float* ua,
    __global const float* va,
    __global float* z
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;    

    int kr = j + (i + 1)*ny;
    int kl = j + (i - 1)*ny;

    int ku = j + 1 + i*ny;  
    int kd = j - 1 + i*ny;  

    if (i == 0) {
		if (j == 0) {
			kr = ny;
			kl = (nx - 1) * ny;

			ku = 1;
			kd = ny - 1;
		}
		else if (j == ny - 1) {
			kr = ny - 1 + ny;
			kl = ny - 1 + (nx - 1) * ny;

			ku = 0;
			kd = ny - 2;
		}
		else {
			kr = j + ny;
			kl = j + (nx - 1) * ny;

			ku = j + 1;
			kd = j - 1;
		}
	}
	else if (i == nx - 1) {
		if (j == 0) {
			kr = 0;
			kl = (nx - 2) * ny;

			ku = 1 + (nx - 1) * ny;
			kd = ny - 1 + (nx - 1) * ny;
		}
		else if (j == ny - 1) {
			kr = ny - 1;
			kl = ny - 1 + (nx - 2) * ny;

			ku = (nx - 1) * ny;
			kd = ny - 2 + (nx - 1) * ny;;
		}
		else {
			kr = j;
			kl = j + (nx - 2) * ny;

			ku = j + 1 + (nx - 1) * ny;
			kd = j - 1 + (nx - 1) * ny;
		}
	}
	else if (j == 0) {
		kr = (i + 1) * ny;
		kl = (i - 1) * ny;

		ku = 1 + i * ny;
		kd = ny - 1 + i * ny;
	}
	else if (j == ny - 1) {
		kr = ny - 1 + (i + 1) * ny;
		kl = ny - 1 + (i - 1) * ny;

		ku = i * ny;
		kd = ny - 2 + i * ny;
	}
    
    z[k] -= dt * (h[kr] * ua[kr] - h[kl] * ua[kl]) / 2 / dx;
	z[k] -= dt * (h[ku] * va[ku] - h[kd] * va[kd]) / 2 / dy;
}

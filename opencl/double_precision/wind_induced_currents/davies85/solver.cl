__kernel void wind_induced_currents_davies85_calc_ua(
    double f, double kb, double g, double rho,
    double dx, double dt,
    int nx, int ny, int nz,
    __global const double* h,
	__global const double* dz,
    __global const double* qx,
    __global const double* ua,
    __global const double* va,
    __global const double* z,
    __global const double* uf,
	__global const double* nu,
    __global double* u1a
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;    
    int idx = nz - 1 + nz*k;

    if (h[k] > 0) {
		int nb = nz + k*(nz + 1);

		double ab = kb*h[k]*dz[nz - 1]/2/nu[nb];
		double ub = (uf[idx] + ua[k])/(1 + ab);

		u1a[k] = ua[k] + f * va[k] * dt - kb * dt * ub / h[k] + qx[k] * dt / rho / h[k];

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

		double gx = (f * va[k] + qx[k] / rho / h[k]) / g;

		double zr = (h[kr] > 0) ? z[kr] : z[k] + dx * gx;
		double zl = (h[kl] > 0) ? z[kl] : z[k] - dx * gx;

		u1a[k] -= g * dt * (zr - zl) / 2 / dx;
	}
}

__kernel void wind_induced_currents_davies85_calc_va(
    double f, double kb, double g, double rho,
    double dy, double dt,
    int nx, int ny, int nz,
    __global const double* h,
	__global const double* dz,
    __global const double* qy,
    __global const double* ua,
    __global const double* va,
    __global const double* z,
    __global const double* vf,
	__global const double* nu,
    __global double* v1a
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;    
    int idx = nz - 1 + nz*k;

    if (h[k] > 0) {
		int nb = nz + k*(nz + 1);

		double ab = kb*h[k]*dz[nz - 1]/2/nu[nb];
		double vb = (vf[idx] + va[k])/(1 + ab);

		v1a[k] = va[k] - f * ua[k] * dt - kb * dt * vb / h[k] + qy[k] * dt / rho / h[k];

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

		double gy = (-f * ua[k] + qy[k] / rho / h[k]) / g;

		double zu = (h[ku] > 0) ? z[ku] : z[k] + dy * gy;
		double zd = (h[kd] > 0) ? z[kd] : z[k] - dy * gy;

		v1a[k] -= g * dt * (zu - zd) / 2 / dy;
	}
}

__kernel void wind_induced_currents_davies85_calc_rhs(
    double f, double kb, double rho, 
    double dx, double dy, double dt,
    int ny, int nz,
    __global const double* nu,
    __global const double* dz,
    __global const double* h,
    __global const double* qx,
    __global const double* qy,
    __global const double* ua,
    __global const double* va,
    __global const double* uf,
    __global const double* vf,
    __global double* ud,
    __global double* vd
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

		int nb = nz + p*(nz + 1);

		double ab = kb*h[p]*dz[nz - 1]/2/nu[nb];

		double ub = (uf[idb] + ua[p])/(1 + ab);
		double vb = (vf[idb] + va[p])/(1 + ab);

		int nid = k + p * (nz + 1);
		int nidu = k + 1 + p * (nz + 1);

		if (k == 0) {
			double um1 = uf[id] + h[p] * qx[p] * dz[k] / rho / nu[nid];
			double vm1 = vf[id] + h[p] * qy[p] * dz[k] / rho / nu[nid];

			ud[id] = uf[id] + f * vf[id] * dt + dt * (
				nu[nidu] * (uf[idu] - uf[id]) / (dz[k + 1] + dz[k]) -
				nu[nid] * (uf[id] - um1) / 2 / dz[k]
			) / h[p] / h[p] / dz[k] + kb * ub * dt / h[p] - qx[p] * dt / rho / h[p]
			+ dt * qx[p] / rho / h[p] / 2 / dz[k];

			vd[id] = vf[id] - f * uf[id] * dt + dt * (
				nu[nidu] * (vf[idu] - vf[id]) / (dz[k + 1] + dz[k]) -
				nu[nid] * (vf[id] - vm1) / 2 / dz[k]
			) / h[p] / h[p] / dz[k] + kb * vb * dt / h[p] - qy[p] * dt / rho / h[p]
			+ dt * qy[p] / rho / h[p] / 2 / dz[k];
		}
		else if (k == nz - 1) {
			double un = (
				uf[id] - kb * h[p] * dz[k] * uf[id] / 2 / nu[nidu]
				- kb * h[p] * dz[k] * ua[p] / nu[nidu]
			) / (1 + kb * h[p] * dz[k] / 2 / nu[nidu]);

			double vn = (
				vf[id] - kb * h[p] * dz[k] * vf[id] / 2 / nu[nidu]
				- kb * h[p] * dz[k] * va[p] / nu[nidu]
			) / (1 + kb * h[p] * dz[k] / 2 / nu[nidu]);

			ud[id] = uf[id] + f * vf[id] * dt + dt * (
				nu[nidu] * (un - uf[id]) / 2 / dz[k] -
				nu[nid] * (uf[id] - uf[idd]) / (dz[k] + dz[k - 1])
			) / h[p] / h[p] / dz[k] + kb * ub * dt / h[p] - qx[p] * dt / rho / h[p]
			- dt * kb * ua[p] / h[p] / 2 / dz[k] / (1 + kb * h[p] * dz[k] / 2 / nu[nidu]);

			vd[id] = vf[id] - f * uf[id] * dt + dt * (
				nu[nidu] * (vn - vf[id]) / 2 / dz[k] -
				nu[nid] * (vf[id] - vf[idd]) / (dz[k] + dz[k - 1])
			) / h[p] / h[p] / dz[k] + kb * vb * dt / h[p] - qy[p] * dt / rho / h[p]
			- dt * kb * va[p] / h[p] / 2 / dz[k] / (1 + kb * h[p] * dz[k] / 2 / nu[nidu]);
		}
		else {
			ud[id] = uf[id] + f * vf[id] * dt + dt * (
				nu[nidu] * (uf[idu] - uf[id]) / (dz[k + 1] + dz[k]) -
				nu[nid] * (uf[id] - uf[idd]) / (dz[k] + dz[k - 1])
			) / h[p] / h[p] / dz[k] + kb * ub * dt / h[p] - qx[p] * dt / rho / h[p];

			vd[id] = vf[id] - f * uf[id] * dt + dt * (
				nu[nidu] * (vf[idu] - vf[id]) / (dz[k + 1] + dz[k]) -
				nu[nid] * (vf[id] - vf[idd]) / (dz[k] + dz[k - 1])
			) / h[p] / h[p] / dz[k] + kb * vb * dt / h[p] - qy[p] * dt / rho / h[p];
		}        
    }    
}

__kernel void wind_induced_currents_davies85_create_tridiagonal_matrix(
    double kb, double dt,
    int ny, int nz,
    __global const double* nu,
    __global const double* dz,
    __global const double* h,
    __global double* al,
    __global double* ac,
    __global double* ar
) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    int k = get_global_id(2);

    int p = j + i*ny;

    if (h[p] > 0) {
		int id = k + p * nz;
		
		int nid = k + p * (nz + 1);
		int nidu = k + 1 + p * (nz + 1);

		if (k == 0) {
			ac[id] = 1 + dt * nu[nidu] / h[p] / h[p] / dz[k] / (dz[k + 1] + dz[k]);
			ar[id] = -dt * nu[nidu] / h[p] / h[p] / dz[k] / (dz[k + 1] + dz[k]);
		}
		else if (k == nz - 1) {
			al[id] = -dt * nu[nid] / h[p] / h[p] / dz[k] / (dz[k] + dz[k - 1]);
			ac[id] = 1 + dt * (
				kb * h[p] * dz[k] / (1 + kb * h[p] * dz[k] / 2 / nu[nidu]) / 2 / dz[k] +
				nu[nid] / (dz[k] + dz[k - 1])
			) / h[p] / h[p] / dz[k];
		}
		else {
			al[id] = -dt * nu[nid] / h[p] / h[p] / dz[k] / (dz[k] + dz[k - 1]);
			ac[id] = 1 + dt * (
				nu[nidu] / (dz[k + 1] + dz[k]) +
				nu[nid] / (dz[k] + dz[k - 1])
			) / h[p] / h[p] / dz[k];
			ar[id] = -dt * nu[nidu] / h[p] / h[p] / dz[k] / (dz[k + 1] + dz[k]);
		}
	}
}

__kernel void wind_induced_currents_davies85_calc_uvf(
    int ny, int nz,
	__global const double* h,
	__global const double* dz,
    __global const double* al,
    __global const double* ac,
    __global const double* ar,
    __global double* uf,
    __global double* vf,
    __global double* ud,
    __global double* vd
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
		
		double mu = 0, mv = 0, s = 0;

		for (int l = 0; l < nz; l++) {
            mu += dz[l] * uf[shift + l];
            mv += dz[l] * vf[shift + l];
            s  += dz[l];
        }

		mu /= s;
        mv /= s;

		for (int l = 0; l < nz; l++) {
            uf[shift + l] -= mu;
            vf[shift + l] -= mv;
        }
    }    
}

__kernel void wind_induced_currents_davies85_calc_residual(
	int ny, int nz,
	__global const double* h,
    __global const double* al,
    __global const double* ac,
    __global const double* ar,
    __global const double* uf,
    __global const double* vf,
    __global double* ud,
    __global double* vd
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
    double dx, double dy, double dt,
    int nx, int ny, int nz,
    __global const double* h,
    __global const double* ua,
    __global const double* va,
    __global double* z
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;    

    if (h[k] > 0) {
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

		double qr = (h[kr] > 0) ? h[kr] * ua[kr] : -h[k] * ua[k];
        double ql = (h[kl] > 0) ? h[kl] * ua[kl] : -h[k] * ua[k];
        double qu = (h[ku] > 0) ? h[ku] * va[ku] : -h[k] * va[k];
        double qd = (h[kd] > 0) ? h[kd] * va[kd] : -h[k] * va[k];

        z[k] -= dt * (qr - ql) / 2 / dx;
        z[k] -= dt * (qu - qd) / 2 / dy;
    }
}

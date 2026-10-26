__kernel void wind_induced_currents_davies85_calc_ua(
    float f, float kb, float g, float rho,
    float dx, float dt,
    int nx, int ny, int nz,
    __global const float* h,
	__global const float* dz,
    __global const float* qx,
    __global const float* ua,
    __global const float* va,
    __global const float* z,
    __global const float* uf,
	__global const float* nu,
    __global float* u1a
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;    
    int idx = nz - 1 + nz*k;

    if (h[k] > 0) {
		int nb = nz + k*(nz + 1);

		float ab = kb*h[k]*dz[nz - 1]/2/nu[nb];
		float ub = (uf[idx] + ua[k])/(1 + ab);

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

		float gx = (f * va[k] + qx[k] / rho / h[k]) / g;

		float zr = (h[kr] > 0) ? z[kr] : z[k] + dx * gx;
		float zl = (h[kl] > 0) ? z[kl] : z[k] - dx * gx;

		u1a[k] -= g * dt * (zr - zl) / 2 / dx;
	}
}

__kernel void wind_induced_currents_davies85_calc_va(
    float f, float kb, float g, float rho,
    float dy, float dt,
    int nx, int ny, int nz,
    __global const float* h,
	__global const float* dz,
    __global const float* qy,
    __global const float* ua,
    __global const float* va,
    __global const float* z,
    __global const float* vf,
	__global const float* nu,
    __global float* v1a
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;    
    int idx = nz - 1 + nz*k;

    if (h[k] > 0) {
		int nb = nz + k*(nz + 1);

		float ab = kb*h[k]*dz[nz - 1]/2/nu[nb];
		float vb = (vf[idx] + va[k])/(1 + ab);

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

		float gy = (-f * ua[k] + qy[k] / rho / h[k]) / g;

		float zu = (h[ku] > 0) ? z[ku] : z[k] + dy * gy;
		float zd = (h[kd] > 0) ? z[kd] : z[k] - dy * gy;

		v1a[k] -= g * dt * (zu - zd) / 2 / dy;
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

__kernel void wind_induced_currents_davies85_calc_rhs_u(
    float f, float kb, float rho, float dt,
    int ny, int nz,
    __global const float* nu,
    __global const float* dz,
    __global const float* h,
    __global const float* qx,
    __global const float* ua,
    __global const float* uf,
    __global const float* vf,
    __global float* ud
) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    int k = get_global_id(2);

    int p = j + i*ny;

    if (h[p] > 0) {
        int id  = k + p*nz;
        int idu = k + 1 + p*nz;
        int idd = k - 1 + p*nz;
        int idb = nz - 1 + p*nz;

        int nid  = k + p*(nz + 1);
        int nidu = k + 1 + p*(nz + 1);
        int nb   = nz + p*(nz + 1);

        float ab = kb*h[p]*dz[nz - 1]/2/nu[nb];
        float ub = (uf[idb] + ua[p])/(1 + ab);

        float visc;

        if (k == 0) {
            float um1 = uf[id] + h[p]*qx[p]*dz[k]/rho/nu[nid];
            visc = nu[nidu]*(uf[idu] - uf[id])/(dz[k + 1] + dz[k]) - nu[nid]*(uf[id] - um1)/2/dz[k];
        } else if (k == nz - 1) {
            float un = (uf[id] - ab*uf[id] - 2*ab*ua[p])/(1 + ab);
            visc = nu[nidu]*(un - uf[id])/2/dz[k] - nu[nid]*(uf[id] - uf[idd])/(dz[k] + dz[k - 1]);
        } else {
            visc = nu[nidu]*(uf[idu] - uf[id])/(dz[k + 1] + dz[k]) - nu[nid]*(uf[id] - uf[idd])/(dz[k] + dz[k - 1]);
        }

        ud[id] = uf[id] + f*vf[id]*dt + dt*visc/h[p]/h[p]/dz[k]
               + kb*ub*dt/h[p] - qx[p]*dt/rho/h[p];

        if (k == 0)      ud[id] += dt*qx[p]/rho/h[p]/2/dz[k];
        if (k == nz - 1) ud[id] -= dt*kb*ua[p]/h[p]/2/dz[k]/(1 + ab);
    }
}

__kernel void wind_induced_currents_davies85_calc_rhs_v(
    float f, float kb, float rho, float dt,
    int ny, int nz,
    __global const float* nu,
    __global const float* dz,
    __global const float* h,
    __global const float* qy,
    __global const float* va,
    __global const float* uf,
    __global const float* vf,    
    __global float* vd
) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    int k = get_global_id(2);

    int p = j + i*ny;

    if (h[p] > 0) {
        int id  = k + p*nz;
        int idu = k + 1 + p*nz;
        int idd = k - 1 + p*nz;
        int idb = nz - 1 + p*nz;

        int nid  = k + p*(nz + 1);
        int nidu = k + 1 + p*(nz + 1);
        int nb   = nz + p*(nz + 1);

        float ab = kb*h[p]*dz[nz - 1]/2/nu[nb];
        float vb = (vf[idb] + va[p])/(1 + ab);

        float visc;

        if (k == 0) {
            float vm1 = vf[id] + h[p]*qy[p]*dz[k]/rho/nu[nid];
            visc = nu[nidu]*(vf[idu] - vf[id])/(dz[k + 1] + dz[k]) - nu[nid]*(vf[id] - vm1)/2/dz[k];
        } else if (k == nz - 1) {
            float vn = (vf[id] - ab*vf[id] - 2*ab*va[p])/(1 + ab);
            visc = nu[nidu]*(vn - vf[id])/2/dz[k] - nu[nid]*(vf[id] - vf[idd])/(dz[k] + dz[k - 1]);
        } else {
            visc = nu[nidu]*(vf[idu] - vf[id])/(dz[k + 1] + dz[k]) - nu[nid]*(vf[id] - vf[idd])/(dz[k] + dz[k - 1]);
        }

        vd[id] = vf[id] - f*uf[id]*dt + dt*visc/h[p]/h[p]/dz[k]
               + kb*vb*dt/h[p] - qy[p]*dt/rho/h[p];

        if (k == 0)      vd[id] += dt*qy[p]/rho/h[p]/2/dz[k];
        if (k == nz - 1) vd[id] -= dt*kb*va[p]/h[p]/2/dz[k]/(1 + ab);
    }
}

__kernel void wind_induced_currents_davies85_calc_column(
    int ny, int nz,
    __global const float* h,
    __global const float* dz,
    __global const float* al,
    __global const float* ac,
    __global const float* ar,
    __global float* d,
    __global float* x
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;

    if (h[k] > 0) {
        int shift = nz*k;

        slae_direct_tridiagonal_calc(nz, al + shift, ac + shift, ar + shift, d + shift, x + shift);

        float m = 0, s = 0;

        for (int l = 0; l < nz; l++) { 
            m += dz[l]*x[shift + l]; 
            s += dz[l]; 
        }

        m /= s;

        for (int l = 0; l < nz; l++) x[shift + l] -= m;
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
    
		float qr = (h[kr] > 0) ? h[kr] * ua[kr] : -h[k] * ua[k];
        float ql = (h[kl] > 0) ? h[kl] * ua[kl] : -h[k] * ua[k];
        float qu = (h[ku] > 0) ? h[ku] * va[ku] : -h[k] * va[k];
        float qd = (h[kd] > 0) ? h[kd] * va[kd] : -h[k] * va[k];

        z[k] -= dt * (qr - ql) / 2 / dx;
        z[k] -= dt * (qu - qd) / 2 / dy;    
    }
}

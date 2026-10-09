#include <iostream>

#include <cublas_v2.h>

#include <parallel/slae/direct/tridiagonal.h>

#include "kernel.h"

using namespace Parallel::SLAE::Direct;
using namespace Parallel::WindInducedCurrents::Davies85;

Kernel::Kernel(
	size_t nx, size_t ny,
	double dx, double dy,
	double rho, double f, double kb, double g,
	vector<double>& dz, vector<double>& h,
	vector<double>& ua, vector<double>& va,
	vector<double>& z,
	vector<double>& uf, vector<double>& vf,
	vector<double>& u10, vector<double>& v10,
	vector<double>& qx, vector<double>& qy,
	vector<double>& nu
): nx(nx), ny(ny),
	dx(dx), dy(dy),
	rho(rho), f(f), kb(kb), g(g),
	dz(dz), h(h),
	ua(ua), va(va),
	z(z),
	uf(uf), vf(vf),
	u10(u10), v10(v10),
	qx(qx), qy(qy),
	nu(nu)
{
	nz = dz.size();

	threadsPerBlock = dim3(8, 8, 8);

	blocksPerGrid2D = dim3(
		(nx + threadsPerBlock.x - 1) / threadsPerBlock.x,
		(ny + threadsPerBlock.y - 1) / threadsPerBlock.y
	);

	blocksPerGrid3D = dim3(
		(nx + threadsPerBlock.x - 1) / threadsPerBlock.x,
		(ny + threadsPerBlock.y - 1) / threadsPerBlock.y,
		(nz + threadsPerBlock.z - 1) / threadsPerBlock.z
	);

	cudaMalloc(&d_dz, nz * sizeof(double));
	cudaMemcpy(d_dz, dz.data(), nz * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_h, nx * ny * sizeof(double));
	cudaMemcpy(d_h, h.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_ua, nx * ny * sizeof(double));
	cudaMemcpy(d_ua, ua.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_u1a, nx * ny * sizeof(double));

	cudaMalloc(&d_va, nx * ny * sizeof(double));
	cudaMemcpy(d_va, va.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_v1a, nx * ny * sizeof(double));

	cudaMalloc(&d_z, nx * ny * sizeof(double));
	cudaMemcpy(d_z, z.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_zt, nx * ny * sizeof(double));
	cudaMemcpy(d_zt, z.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_uf, nx * ny * nz * sizeof(double));
	cudaMemcpy(d_uf, uf.data(), nx * ny * nz * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_uft, nx * ny * nz * sizeof(double));
	cudaMemcpy(d_uft, uf.data(), nx * ny * nz * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_ud, nx * ny * nz * sizeof(double));
	cudaMalloc(&d_udt, nx * ny * nz * sizeof(double));

	cudaMalloc(&d_vf, nx * ny * nz * sizeof(double));
	cudaMemcpy(d_vf, vf.data(), nx * ny * nz * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_vft, nx * ny * nz * sizeof(double));
	cudaMemcpy(d_vft, vf.data(), nx * ny * nz * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_vd, nx * ny * nz * sizeof(double));
	cudaMalloc(&d_vdt, nx * ny * nz * sizeof(double));

	cudaMalloc(&d_al, nx * ny * nz * sizeof(double));
	cudaMalloc(&d_ac, nx * ny * nz * sizeof(double));
	cudaMalloc(&d_ar, nx * ny * nz * sizeof(double));

	cudaMalloc(&d_u10, nx * ny * sizeof(double));
	cudaMemcpy(d_u10, u10.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_v10, nx * ny * sizeof(double));
	cudaMemcpy(d_v10, v10.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_qx, nx * ny * sizeof(double));
	cudaMemcpy(d_qx, qx.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_qy, nx * ny * sizeof(double));
	cudaMemcpy(d_qy, qy.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_nu, nx * ny * (nz + 1) * sizeof(double));
	cudaMemcpy(d_nu, nu.data(), nx * ny * (nz + 1) * sizeof(double), cudaMemcpyHostToDevice);

	double alpha = -1.0;

	cudaMalloc(&d_alpha, sizeof(double));
	cudaMemcpy(d_alpha, &alpha, sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_res, sizeof(double));
	cudaMalloc(&d_res2, sizeof(double));

	cudaMalloc(&d_res_idx, sizeof(int));
}

__global__ void __extract_element__(const double* src, const int* index, double* dst) {
	if (*index > 0) {
		size_t i = *index - 1;

		*dst = src[i];
	}
}

__global__ void __calcUA__(
	const double rho, const double f, const double kb, const double g,
	const double dx, const double dt,
	const size_t nx, const size_t ny, const size_t nz,
	const double* h, 
	const double* qx,
	const double* ua,
	const double* va,
	const double* z,
	const double* uf,
	double* u1a
) {
	size_t i = blockIdx.x * blockDim.x + threadIdx.x;
	size_t j = blockIdx.y * blockDim.y + threadIdx.y;

	if (i < nx && j < ny) {
		size_t k = j + i * ny;
		size_t idx = nz - 1 + nz * k;

		if (h[k] > 0) {
			u1a[k] = ua[k] + f * va[k] * dt - kb * dt * (uf[idx] + ua[k]) / h[k] + qx[k] * dt / rho / h[k];

			if (i > 0 && i < nx - 1) {
				size_t kr = j + (i + 1) * ny;
				size_t kl = j + (i - 1) * ny;


				u1a[k] -= g * dt * (z[kr] - z[kl]) / 2 / dx;
			}
		}
	}
}

void Kernel::calcUA(double dt) {
	__calcUA__ << < blocksPerGrid2D, threadsPerBlock >> > (rho, f, kb, g, dx, dt, nx, ny, nz, d_h, d_qx, d_ua, d_va, d_z, d_uf, d_u1a);

	cudaDeviceSynchronize();

	cudaError_t err = cudaGetLastError();

	if (err != cudaSuccess) {
		cerr << "CUDA Error: " << cudaGetErrorString(err) << endl;
	}
}

__global__ void __calcVA__(
	const double rho, const double f, const double kb, const double g,
	const double dy, const double dt,
	const size_t nx, const size_t ny, const size_t nz,
	const double* h,
	const double* qy,
	const double* ua,
	const double* va,
	const double* z,
	const double* vf,
	double* v1a
) {
	size_t i = blockIdx.x * blockDim.x + threadIdx.x;
	size_t j = blockIdx.y * blockDim.y + threadIdx.y;

	if (i < nx && j < ny) {
		size_t k = j + i * ny;
		size_t idx = nz - 1 + nz * k;

		if (h[k] > 0) {
			v1a[k] = va[k] - f * ua[k] * dt - kb * dt * (vf[idx] + va[k]) / h[k] + qy[k] * dt / rho / h[k];

			if (j > 0 && j < ny - 1) {
				size_t ku = j + 1 + i * ny;
				size_t kd = j - 1 + i * ny;

				v1a[k] -= g * dt * (z[ku] - z[kd]) / 2 / dy;
			}
		}
	}
}

void Kernel::calcVA(double dt) {
	__calcVA__ << < blocksPerGrid2D, threadsPerBlock >> > (rho, f, kb, g, dy, dt, nx, ny, nz, d_h, d_qy, d_ua, d_va, d_z, d_vf, d_v1a);

	cudaDeviceSynchronize();

	cudaError_t err = cudaGetLastError();

	if (err != cudaSuccess) {
		cerr << "CUDA Error: " << cudaGetErrorString(err) << endl;
	}
}

void Kernel::updateUVA() {
	cudaMemcpy(ua.data(), d_ua, nx * ny * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(va.data(), d_va, nx * ny * sizeof(double), cudaMemcpyDeviceToHost);
}

void Kernel::swapUA() {
	swap(d_ua, d_u1a);
}

void Kernel::swapVA() {	
	swap(d_va, d_v1a);
}

__global__ void __calcRHS__(
	const double rho, const double f, const double kb,
	const double dx, const double dy, const double dt,
	const size_t nx, const size_t ny, const size_t nz,
	const double* dz, const double* h,	
	const double* qx, const double* qy,
	const double* ua, const double* u1a,
	const double* va, const double* v1a,
	const double* uf, const double* vf,
	const double* nu,
	double* ud, double* vd
) {
	size_t i = blockIdx.x * blockDim.x + threadIdx.x;
	size_t j = blockIdx.y * blockDim.y + threadIdx.y;
	size_t k = blockIdx.z * blockDim.z + threadIdx.z;

	if (i < nx && j < ny && k < nz) {
		size_t p = j + i * ny;

		size_t id = k + p * nz;

		if (h[p] > 0) {
			size_t idu = k + 1 + p * nz;
			size_t idd = k - 1 + p * nz;

			size_t idb = nz - 1 + p * nz;

			if (k == 0) {
				double um1 = uf[id] + h[p] * qx[p] * dz[k] / rho / nu[id];
				double vm1 = vf[id] + h[p] * qy[p] * dz[k] / rho / nu[id];

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
				double un = (
					uf[id] - kb * h[p] * dz[k] * uf[id] / 2 / nu[idu]
					- kb * h[p] * dz[k] * ua[p] / nu[idu]
				) / (1 + kb * h[p] * dz[k] / 2 / nu[idu]);

				double vn = (
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
				double u = uf[idb] + ua[p];
				double v = vf[idb] + va[p];

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
}

void Kernel::calcRHS(double dt) {
	__calcRHS__ << < blocksPerGrid3D, threadsPerBlock >> > (rho, f, kb, dx, dy, dt, nx, ny, nz, d_dz, d_h, d_qx, d_qy, d_ua, d_u1a, d_va, d_v1a, d_uf, d_vf, d_nu, d_ud, d_vd);

	cudaDeviceSynchronize();

	cudaError_t err = cudaGetLastError();

	if (err != cudaSuccess) {
		cerr << "CUDA Error: " << cudaGetErrorString(err) << endl;
	}

	cudaMemcpy(d_udt, d_ud, nx * ny * nz * sizeof(double), cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_vdt, d_vd, nx * ny * nz * sizeof(double), cudaMemcpyDeviceToDevice);
}

__global__ void __createTridiagonalMatrix__(
	const double kb, const double dt,
	const size_t nx, const size_t ny, const size_t nz,
	const double* dz, const double* h,
	const double* nu,
	double* al, double* ac, double* ar
) {
	size_t i = blockIdx.x * blockDim.x + threadIdx.x;
	size_t j = blockIdx.y * blockDim.y + threadIdx.y;
	size_t k = blockIdx.z * blockDim.z + threadIdx.z;

	if (i < nx && j < ny && k < nz) {
		size_t p = j + i * ny;

		size_t id = k + p * nz;
		size_t idu = k + 1 + p * nz;

		if (h[p] > 0) {
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
}

void Kernel::createTridiagonalMatrix(double dt) {
	__createTridiagonalMatrix__ << < blocksPerGrid3D, threadsPerBlock >> > (kb, dt, nx, ny, nz, d_dz, d_h, d_nu, d_al, d_ac, d_ar);

	cudaDeviceSynchronize();

	cudaError_t err = cudaGetLastError();

	if (err != cudaSuccess) {
		cerr << "CUDA Error: " << cudaGetErrorString(err) << endl;
	}
}

__global__ void __calcUVF__(
	const size_t nx, const size_t ny, const size_t nz,
	const double* h,
	const double* al, const double* ac, const double* ar,
	double* uf, double* vf,
	double* ud, double* vd
) {
	size_t i = blockIdx.x * blockDim.x + threadIdx.x;
	size_t j = blockIdx.y * blockDim.y + threadIdx.y;

	if (i < nx && j < ny) {
		size_t k = j + i * ny;

		if (h[k] > 0) {
			size_t shift = nz * k;

			Tridiagonal::Matrix matrix(nz, al + shift, ac + shift, ar + shift);

			matrix.solve(uf + shift, ud + shift);
			matrix.solve(vf + shift, vd + shift);
		}
	}
}

void Kernel::calcUVF() {
	__calcUVF__ << < blocksPerGrid2D, threadsPerBlock >> > (nx, ny, nz, d_h, d_al, d_ac, d_ar, d_uf, d_vf, d_ud, d_vd);

	cudaDeviceSynchronize();

	cudaError_t err = cudaGetLastError();

	if (err != cudaSuccess) {
		cerr << "CUDA Error: " << cudaGetErrorString(err) << endl;
	}
}

void Kernel::updateUVF() {
	cudaMemcpy(uf.data(), d_uf, nx * ny * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(vf.data(), d_vf, nx * ny * sizeof(double), cudaMemcpyDeviceToHost);
}

void Kernel::updateUVFTDevice() {
	cudaMemcpy(d_uft, d_uf, nx * ny * nz * sizeof(double), cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_vft, d_vf, nx * ny * nz * sizeof(double), cudaMemcpyDeviceToDevice);
}

double Kernel::maxAbsUFDifference() {
	cublasHandle_t handle;
	cublasCreate(&handle);

	cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);

	cublasDaxpy(handle, nx * ny * nz, d_alpha, d_uf, 1, d_uft, 1);
	cublasIdamax(handle, nx * ny * nz, d_uft, 1, d_res_idx);

	__extract_element__ << < 1, 1 >> > (d_uft, d_res_idx, d_res);

	cudaDeviceSynchronize();

	cudaError_t err = cudaGetLastError();

	if (err != cudaSuccess) {
		cerr << "CUDA Error: " << cudaGetErrorString(err) << endl;
	}

	double res;

	cudaMemcpy(&res, d_res, sizeof(double), cudaMemcpyDeviceToHost);

	cublasDestroy(handle);

	return abs(res);
}

double Kernel::maxAbsVFDifference() {
	cublasHandle_t handle;
	cublasCreate(&handle);

	cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);

	cublasDaxpy(handle, nx * ny * nz, d_alpha, d_vf, 1, d_vft, 1);
	cublasIdamax(handle, nx * ny * nz, d_vft, 1, d_res_idx);

	__extract_element__ << < 1, 1 >> > (d_vft, d_res_idx, d_res);

	cudaDeviceSynchronize();

	cudaError_t err = cudaGetLastError();

	if (err != cudaSuccess) {
		cerr << "CUDA Error: " << cudaGetErrorString(err) << endl;
	}

	double res;

	cudaMemcpy(&res, d_res, sizeof(double), cudaMemcpyDeviceToHost);

	cublasDestroy(handle);

	return abs(res);
}

__global__ void __calcZ__(
	const double dx, const double dy, const double dt,
	const size_t nx, const size_t ny,
	const double* h,
	const double* ua, const double* va,
	double* z
) {
	size_t i = blockIdx.x * blockDim.x + threadIdx.x;
	size_t j = blockIdx.y * blockDim.y + threadIdx.y;

	if (i < nx && j < ny) {
		size_t k = j + i * ny;

		if (h[k] > 0) {
			if (i > 0 && i < nx - 1) {
				size_t kr = j + (i + 1) * ny;
				size_t kl = j + (i - 1) * ny;

				z[k] -= dt * (h[kr] * ua[kr] - h[kl] * ua[kl]) / 2 / dx;
			}

			if (j > 0 && j < ny - 1) {
				size_t ku = j + 1 + i * ny;
				size_t kd = j - 1 + i * ny;

				z[k] -= dt * (h[ku] * va[ku] - h[kd] * va[kd]) / 2 / dy;
			}
		}
	}
}

void Kernel::calcZ(double dt) {
	__calcZ__ << < blocksPerGrid2D, threadsPerBlock >> > (dx, dy, dt, nx, ny, d_h, d_ua, d_va, d_z);

	cudaDeviceSynchronize();

	cudaError_t err = cudaGetLastError();

	if (err != cudaSuccess) {
		cerr << "CUDA Error: " << cudaGetErrorString(err) << endl;
	}
}

void Kernel::updateZ() {
	cudaMemcpy(z.data(), d_z, nx * ny * sizeof(double), cudaMemcpyDeviceToHost);
}

void Kernel::updateZTDevice() {
	cudaMemcpy(d_zt, d_z, nx * ny * sizeof(double), cudaMemcpyDeviceToDevice);
}

double Kernel::maxAbsZDifference() {
	cublasHandle_t handle;
	cublasCreate(&handle);

	cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);

	cublasDaxpy(handle, nx * ny * nz, d_alpha, d_z, 1, d_zt, 1);
	cublasIdamax(handle, nx * ny * nz, d_zt, 1, d_res_idx);

	__extract_element__ << < 1, 1 >> > (d_zt, d_res_idx, d_res);

	cudaDeviceSynchronize();

	cudaError_t err = cudaGetLastError();

	if (err != cudaSuccess) {
		cerr << "CUDA Error: " << cudaGetErrorString(err) << endl;
	}

	double res;

	cudaMemcpy(&res, d_res, sizeof(double), cudaMemcpyDeviceToHost);

	cublasDestroy(handle);

	return abs(res);
}

__global__ void __calcResidual__(
	const size_t nx, const size_t ny, const size_t nz, 
	const double* h,
	const double* al, const double* ac, const double* ar,
	const double* uf, const double* vf,
	double* ud, double* vd
) {
	size_t i = blockIdx.x * blockDim.x + threadIdx.x;
	size_t j = blockIdx.y * blockDim.y + threadIdx.y;

	if (i < nx && j < ny) {
		size_t k = j + i * ny;

		if (h[k] > 0) {
			size_t shift = nz * k;

			Tridiagonal::Matrix matrix(nz, al + shift, ac + shift, ar + shift);

			matrix.residual(uf + shift, ud + shift);
			matrix.residual(vf + shift, vd + shift);
		}
	}
}

double Kernel::maxAbsResidual() {
	__calcResidual__ << < blocksPerGrid2D, threadsPerBlock >> > (nx, ny, nz, d_h, d_al, d_ac, d_ar, d_uf, d_vf, d_udt, d_vdt);

	cublasHandle_t handle;
	cublasCreate(&handle);

	cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);

	cublasIdamax(handle, nx * ny * nz, d_udt, 1, d_res_idx);

	__extract_element__ << < 1, 1 >> > (d_udt, d_res_idx, d_res);

	cublasIdamax(handle, nx * ny * nz, d_vdt, 1, d_res_idx);

	__extract_element__ << < 1, 1 >> > (d_vdt, d_res_idx, d_res2);

	cudaDeviceSynchronize();

	cudaError_t err = cudaGetLastError();

	if (err != cudaSuccess) {
		cerr << "CUDA Error: " << cudaGetErrorString(err) << endl;
	}

	double res, res2;

	cudaMemcpy(&res, d_res, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&res2, d_res2, sizeof(double), cudaMemcpyDeviceToHost);

	cublasDestroy(handle);

	return max(abs(res), abs(res2));
}

void Kernel::updateQUV10Device() {
	cudaMemcpy(d_u10, u10.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_v10, v10.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(d_qx, qx.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_qy, qy.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice);
}

void Kernel::calcNU() {}

void Kernel::updateNU() {}

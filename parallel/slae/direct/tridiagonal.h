#ifndef PARALLEL_SLAE_DIRECT_TRIDIAGONAL_H
#define PARALLEL_SLAE_DIRECT_TRIDIAGONAL_H

#include <cuda_runtime.h>

namespace Parallel {
	namespace SLAE {
		namespace Direct {
			namespace Tridiagonal {
				class Matrix {
				private:
					const size_t n;
					const double* l;
					const double* c;
					const double* r;

				public:
					__device__ Matrix(const size_t n, const double* l, const double* c, const double* r);

					__device__ void solve(double* u, double* d);
					__device__ void residual(const double* u, double* d);
				};
			}
		}
	}
}

#endif
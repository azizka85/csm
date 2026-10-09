#include "tridiagonal.h"

namespace Parallel {
	namespace SLAE {
		namespace Direct {
			namespace Tridiagonal {
				__device__ Matrix::Matrix(const size_t n, const double* l, const double* c, const double* r): n(n) {
					this->l = l;
					this->c = c;
					this->r = r;
				}

				__device__ void Matrix::solve(double* u, double* d) {
					if (n > 0) {
						u[0] = r[0] / c[0];
						d[0] = d[0] / c[0];

						for (size_t i = 1; i < n; i++) {
							auto c1 = c[i] - l[i] * u[i - 1];

							u[i] = r[i] / c1;
							d[i] = (d[i] - l[i] * d[i - 1]) / c1;
						}

						u[n - 1] = d[n - 1];

						for (int i = n - 2; i >= 0; i--) {
							u[i] = d[i] - u[i] * u[i + 1];
						}
					}
				}

				__device__ void Matrix::residual(const double* u, double* d) {
					if (n > 0) {
						if (n == 1) {
							d[0] -= c[0] * u[0];
						}
						else {
							d[0] -= c[0] * u[0] + r[0] * u[1];

							for (size_t i = 1; i < n - 1; i++) {
								d[i] -= c[i] * u[i] + l[i] * u[i - 1] + r[i] * u[i + 1];
							}

							d[n - 1] -= c[n - 1] * u[n - 1] + l[n - 1] * u[n - 2];
						}
					}
				}
			}
		}
	}
}
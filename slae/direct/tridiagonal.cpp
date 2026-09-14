#include <iostream>

#include <utils/checks/vector.h>

#include "tridiagonal.h"

using namespace Utils::Checks;

using namespace SLAE::Direct;

Tridiagonal::UniformMatrix::UniformMatrix(UniformMatrixParams params) {
	this->params = move(params);
}

void Tridiagonal::UniformMatrix::solve(SolutionState state) {
	auto& u = state.u;
	auto& d = state.d;

	const auto l = params.l;
	const auto l1 = params.l1;

	const auto c = params.c;
	const auto c0 = params.c0;
	const auto c1 = params.c1;

	const auto r = params.r;
	const auto r0 = params.r0;

	auto n = d.size();

	Vector::checkU(n, u);

    if (n > 0) {
        u[0] = r0 / c0;
        d[0] = d[0] / c0;

        if (n > 1) {
            for (size_t i = 1; i < n - 1; i++) {
                auto ct = c - l * u[i - 1];

                u[i] = r / ct;
                d[i] = (d[i] - l * d[i - 1]) / ct;
            }

            auto ct = c1 - l1 * u[n - 2];

            d[n - 1] = (d[n - 1] - l1 * d[n - 2]) / ct;
        }

        u[n - 1] = d[n - 1];

        for (int i = n - 2; i >= 0; i--) {
            u[i] = d[i] - u[i] * u[i + 1];
        }
    }
}

void Tridiagonal::UniformMatrix::residual(ResidualState state) {
	auto& u = state.u;
	auto& d = state.d;

	const auto l = params.l;
    const auto l1 = params.l1;

    const auto c = params.c;
    const auto c0 = params.c0;
    const auto c1 = params.c1;

    const auto r = params.r;
    const auto r0 = params.r0;

	auto n = d.size();

	Vector::checkU(n, u);

    if (n > 0) {
        if (n == 1) {
            d[0] -= c0 * u[0];
        }
        else {
            d[0] -= c0 * u[0] + r0 * u[1];

            for (size_t i = 1; i < n - 1; i++) {
                d[i] -= c * u[i] + l * u[i - 1] + r * u[i + 1];
            }

            d[n - 1] -= c1 * u[n - 1] + l1 * u[n - 2];
        }
    }
}

Tridiagonal::Matrix::Matrix(MatrixParams params) {
    Vector::check(params.l, params.c, params.r);

    this->params = move(params);
}

void Tridiagonal::Matrix::solve(SolutionState state) {
    auto& u = state.u;
    auto& d = state.d;

    auto l = params.l;
    auto c = params.c;
    auto r = params.r;

    auto n = c.size();

    Vector::checkU(n, u);
    Vector::checkU(n, d);

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

void Tridiagonal::Matrix::residual(ResidualState state) {
    auto& u = state.u;
    auto& d = state.d;

    auto l = params.l;
    auto c = params.c;
    auto r = params.r;

    auto n = c.size();

    Vector::checkU(n, u);
    Vector::checkU(n, d);

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

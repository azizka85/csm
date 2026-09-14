#define BOOST_TEST_MODULE SLAE_Tridiagonal_Tests

#include <format>

#include <cmath>

#include <stdexcept>

#include <boost/test/included/unit_test.hpp>

#include <slae/direct/tridiagonal.h>

using namespace SLAE::Direct;

BOOST_AUTO_TEST_CASE(test_tridiagonal_matrix) {    
    double epsilon = 1e-13;

    double h = 1;
    double dy = 0.002;

    size_t ny = ceil(h / dy) + 1;

    Tridiagonal::UniformMatrix matrix(Tridiagonal::UniformMatrixParams {
        .l = 1,
        .l1 = 0,

        .c = -2,
        .c0 = 1,
        .c1 = 1,

        .r = 1,
        .r0 = 0
    });

    vector<double> u(ny);
    span<double> us(u);

    vector<double> d(ny-1);
    span<double> ds(d);

    BOOST_CHECK_EXCEPTION(
        matrix.solve(Tridiagonal::SolutionState{
           .u = us,
           .d = ds
        }),
        runtime_error,
        [ny](const runtime_error& e) {
            return string(e.what()) == format(
                "The lengths of the u should be {}, but now the length of the u is {}",
                ny - 1, ny
            );
        }
    );    

    BOOST_CHECK_EXCEPTION(
        matrix.residual(Tridiagonal::ResidualState{
           .u = us,
           .d = ds
        }),
        runtime_error,
        [ny](const runtime_error& e) {
            return string(e.what()) == format(
                "The lengths of the u should be {}, but now the length of the u is {}",
                ny - 1, ny
            );
        }
    );

    vector<double> l(ny - 1);
    vector<double> c(ny);
    vector<double> r(ny - 1);

    BOOST_CHECK_EXCEPTION(
        Tridiagonal::Matrix(Tridiagonal::MatrixParams {
            .l = l,
            .c = c,
            .r = r
        }),
        runtime_error,
        [ny](const runtime_error& e) {
            return string(e.what()) == format(
                "The lengths of the l and c should be equal, but now the length of the l is {} and c is {}",
                ny - 1, ny
            );
        }
    );

    l = vector<double>(ny);

    BOOST_CHECK_EXCEPTION(
        Tridiagonal::Matrix(Tridiagonal::MatrixParams{
            .l = l,
            .c = c,
            .r = r
        }),
        runtime_error,
        [ny](const runtime_error& e) {
            return string(e.what()) == format(
                "The lengths of the r and c should be equal, but now the length of the r is {} and c is {}",
                ny - 1, ny
            );
        }
    );
}

BOOST_AUTO_TEST_CASE(test_laplace_equation_with_dirichlet_bc) {
    double diff_epsilon = 1e-13;
    double res_epsilon = 1e-15;

    double h = 1;
    double dy = 0.002;

    size_t ny = ceil(h / dy) + 1;

    Tridiagonal::UniformMatrix matrix(Tridiagonal::UniformMatrixParams{
        .l = 1,
        .l1 = 0,

        .c = -2,
        .c0 = 1,
        .c1 = 1,

        .r = 1,
        .r0 = 0
    });

    vector<double> u(ny);
    span<double> us(u);

    vector<double> d(ny);

    d[0] = 0;
    d[ny - 1] = 1;

    vector<double> d1 = d;

    span<double> ds(d);

    matrix.solve(Tridiagonal::SolutionState {
        .u = us,
        .d = ds
    });

    for (size_t i = 0; i < ny; i++) {
        auto diff = abs(u[i] - i * dy);

        BOOST_CHECK_LT(diff, diff_epsilon);
    }

    ds = span<double>(d1);

    matrix.residual(Tridiagonal::ResidualState{
        .u = u,
        .d = ds
    });

    for (size_t i = 0; i < ny; i++) {
        auto res = abs(d1[i]);

        BOOST_CHECK_LT(res, res_epsilon);
    }

    vector<double> l(ny, 1);

    l[0] = 0;
    l[ny - 1] = 0;

    vector<double> c(ny, -2);

    c[0] = 1;
    c[ny - 1] = 1;

    vector<double> r(ny, 1);

    r[0] = 0;
    r[ny - 1] = 0;

    Tridiagonal::Matrix matrix2(Tridiagonal::MatrixParams {
        .l = l,
        .c = c,
        .r = r
    });

    d = vector<double>(ny);

    d[0] = 0;
    d[ny - 1] = 1;

    d1 = d;

    ds = span<double>(d);

    matrix2.solve(Tridiagonal::SolutionState {
        .u = us,
        .d = ds
    });

    for (size_t i = 0; i < ny; i++) {
        auto diff = abs(u[i] - i * dy);

        BOOST_CHECK_LT(diff, diff_epsilon);
    }

    ds = span<double>(d1);

    matrix2.residual(Tridiagonal::ResidualState{
        .u = u,
        .d = ds
    });

    for (size_t i = 0; i < ny; i++) {
        auto res = abs(d1[i]);

        BOOST_CHECK_LT(res, res_epsilon);
    }
}

BOOST_AUTO_TEST_CASE(test_laplace_equation_with_neumann_bc) {
    double diff_epsilon = 1e-12;
    double res_epsilon = 1e-15;

    double h = 1;
    double dy = 0.002;

    size_t ny = ceil(h / dy) + 1;

    Tridiagonal::UniformMatrix matrix(Tridiagonal::UniformMatrixParams{
        .l = 1,
        .l1 = 2,

        .c = -2,
        .c0 = 1,
        .c1 = -2,

        .r = 1,
        .r0 = 0
    });

    vector<double> u(ny);
    span<double> us(u);

    vector<double> d(ny);

    d[0] = 1;
    d[ny - 1] = 0;

    vector<double> d1 = d;

    span<double> ds(d);

    matrix.solve(Tridiagonal::SolutionState{
        .u = us,
        .d = ds
    });

    for (size_t i = 0; i < ny; i++) {
        auto diff = abs(u[i] - 1);

        BOOST_CHECK_LT(diff, diff_epsilon);
    }

    ds = span<double>(d1);

    matrix.residual(Tridiagonal::ResidualState{
        .u = u,
        .d = ds
    });

    for (size_t i = 0; i < ny; i++) {
        auto res = abs(d1[i]);

        BOOST_CHECK_LT(res, res_epsilon);
    }

    vector<double> l(ny, 1);

    l[0] = 0;
    l[ny - 1] = 2;

    vector<double> c(ny, -2);

    c[0] = 1;
    c[ny - 1] = -2;

    vector<double> r(ny, 1);

    r[0] = 0;
    r[ny - 1] = 0;

    Tridiagonal::Matrix matrix2(Tridiagonal::MatrixParams{
        .l = l,
        .c = c,
        .r = r
    });

    d = vector<double>(ny);

    d[0] = 1;
    d[ny - 1] = 0;

    d1 = d;

    ds = span<double>(d);

    matrix2.solve(Tridiagonal::SolutionState{
        .u = us,
        .d = ds
    });

    for (size_t i = 0; i < ny; i++) {
        auto diff = abs(u[i] - 1);

        BOOST_CHECK_LT(diff, diff_epsilon);
    }

    ds = span<double>(d1);

    matrix2.residual(Tridiagonal::ResidualState{
        .u = u,
        .d = ds
    });

    for (size_t i = 0; i < ny; i++) {
        auto res = abs(d1[i]);

        BOOST_CHECK_LT(res, res_epsilon);
    }
}
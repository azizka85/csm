#include <format>

#include <calc/grid.h>

using namespace Calc;

BOOST_AUTO_TEST_CASE(test_triple_point_vertical_step_size_params) {
    Grid::TriplePointVerticalStepSizeParams params {
        .zm = 0.3,
        .dzMin = -0.1,
        .dzMax = -0.1
    };

    BOOST_CHECK_EXCEPTION(
        Grid::TriplePointVerticalStepSize dz(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("zm should be in the range [0.5, 1], but it is {}", params.zm);
        }
    );

    params.zm = 1.2;

    BOOST_CHECK_EXCEPTION(
        Grid::TriplePointVerticalStepSize dz(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("zm should be in the range [0.5, 1], but it is {}", params.zm);
        }
    );

    params.zm = 0.7;

    BOOST_CHECK_EXCEPTION(
        Grid::TriplePointVerticalStepSize dz(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("dz should be > 0, but it is {}", params.dzMin);
        }
    );

    params.dzMin = 0.1;

    BOOST_CHECK_EXCEPTION(
        Grid::TriplePointVerticalStepSize dz(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("dz should be > 0, but it is {}", params.dzMax);
        }
    );

    params.dzMax = 0.1;

    BOOST_CHECK_EXCEPTION(
        Grid::TriplePointVerticalStepSize dz(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("dzMin should be < dzMax, but they are {} and {}", params.dzMin, params.dzMax);
        }
    );
}

BOOST_AUTO_TEST_CASE(test_triple_point_vertical_step_size_calc) {
    Grid::TriplePointVerticalStepSizeParams params{
        .zm = 0.51,
        .dzMin = 0.002,
        .dzMax = 0.1
    };

    Grid::TriplePointVerticalStepSize dz(params);

    BOOST_CHECK_EQUAL(
        dz.dirName("tpdz"), 
        format(
            "tpdz, dz={}-{}, z={}",
            params.dzMin, params.dzMax, params.zm
        )
    );

    BOOST_CHECK_EQUAL(dz.calculate(0), params.dzMin);

    BOOST_CHECK_EQUAL(dz.calculate(params.zm), params.dzMin + params.dzMax);

    BOOST_CHECK_GT(dz.calculate(1), params.dzMin);
}

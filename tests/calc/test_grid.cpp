#include <format>

#include <doctest/doctest.h>

#include <calc/grid.h>

using namespace Calc;

TEST_CASE("Test triple point vertical step size params") {
    Grid::TriplePointVerticalStepSizeParams params {
        .zm = 0.3,
        .dzMin = -0.1,
        .dzMax = -0.1
    };

    CHECK_THROWS_WITH_AS(
        Grid::TriplePointVerticalStepSize dz(params),
        format("zm should be in the range [0.5, 1], but it is {}", params.zm).c_str(),
        runtime_error
    );

    params.zm = 1.2;

    CHECK_THROWS_WITH_AS(
        Grid::TriplePointVerticalStepSize dz(params),
        format("zm should be in the range [0.5, 1], but it is {}", params.zm).c_str(),
        runtime_error
    );

    params.zm = 0.7;

    CHECK_THROWS_WITH_AS(
        Grid::TriplePointVerticalStepSize dz(params),
        format("dz should be > 0, but it is {}", params.dzMin).c_str(),
        runtime_error
    );

    params.dzMin = 0.1;

    CHECK_THROWS_WITH_AS(
        Grid::TriplePointVerticalStepSize dz(params),
        format("dz should be > 0, but it is {}", params.dzMax).c_str(),
        runtime_error
    );

    params.dzMax = 0.1;

    CHECK_THROWS_WITH_AS(
        Grid::TriplePointVerticalStepSize dz(params),
        format("dzMin should be < dzMax, but they are {} and {}", params.dzMin, params.dzMax).c_str(),
        runtime_error
    );
}

TEST_CASE("Test triple point vertical step size calc") {
    Grid::TriplePointVerticalStepSizeParams params{
        .zm = 0.51,
        .dzMin = 0.002,
        .dzMax = 0.1
    };

    Grid::TriplePointVerticalStepSize dz(params);

    CHECK_EQ(
        dz.dirName("tpdz"), 
        format(
            "tpdz, dz={}-{}, z={}",
            params.dzMin, params.dzMax, params.zm
        )
    );

    CHECK_EQ(dz.calculate(0), params.dzMin);

    CHECK_EQ(dz.calculate(params.zm), params.dzMin + params.dzMax);

    CHECK_GT(dz.calculate(1), params.dzMin);
}

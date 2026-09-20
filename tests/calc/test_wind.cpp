#include <format>

#include <doctest/doctest.h>

#include <calc/wind.h>

using namespace std;

using namespace Calc;

TEST_CASE("Test default wind stress params") {
    Wind::DefaultStressParams params{
        .rhoAir = 0,
        .cd = -1
    };

    CHECK_THROWS_WITH_AS(
        Wind::DefaultStress wind(params),
        format("rho should be > 0, but it is {}", params.rhoAir).c_str(),
        runtime_error
    );

    params.rhoAir = 1.225;

    CHECK_THROWS_WITH_AS(
        Wind::DefaultStress wind(params),
        format("cd should be > 0, but it is {}", params.cd).c_str(),
        runtime_error
    );
}

TEST_CASE("Test default wind stress calc") {
    double epsilon = 1e-15;

    Wind::DefaultStressParams params{
        .rhoAir = 1.225,
        .cd = 2.5e-3
    };

    Wind::DefaultStress wind(params);

    CHECK_EQ(
        wind.dirName("wind"),
        format("wind, rho_air={}, cd={}", params.rhoAir, params.cd)
    );

    auto stress1 = wind.calculate(Wind::SpeedVector {
        .u10 = 22,
        .v10 = 22
    });

    auto stressMagnitude1 = sqrt(stress1.qx*stress1.qx + stress1.qy*stress1.qy);

    CHECK_LT(abs(stressMagnitude1 - 2.9645), epsilon);

    auto stress2 = wind.calculate(Wind::SpeedVector {
        .u10 = 22,
        .v10 = 0
    });

    CHECK_EQ(stress2.qx, 1.48225);
    CHECK_EQ(stress2.qy, 0);
}
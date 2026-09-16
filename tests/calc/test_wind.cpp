#include <calc/wind.h>

using namespace Calc;

BOOST_AUTO_TEST_CASE(test_default_wind_stress_params) {
    Wind::DefaultStressParams params{
        .rhoAir = 0,
        .cd = -1
    };

    BOOST_CHECK_EXCEPTION(
        Wind::DefaultStress wind(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("rho should be > 0, but it is {}", params.rhoAir);
        }
    );

    params.rhoAir = 1.225;

    BOOST_CHECK_EXCEPTION(
        Wind::DefaultStress wind(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("cd should be > 0, but it is {}", params.cd);
        }
    );
}

BOOST_AUTO_TEST_CASE(test_default_wind_stress_calc) {
    double epsilon = 1e-15;

    Wind::DefaultStressParams params{
        .rhoAir = 1.225,
        .cd = 2.5e-3
    };

    Wind::DefaultStress wind(params);

    BOOST_CHECK_EQUAL(
        wind.dirName("wind"),
        format("wind, rho_air={}, cd={}", params.rhoAir, params.cd)
    );

    auto stress1 = wind.calculate(Wind::SpeedVector {
        .u10 = 22,
        .v10 = 22
    });

    auto stressMagnitude1 = sqrt(stress1.qx*stress1.qx + stress1.qy*stress1.qy);

    BOOST_CHECK_LT(abs(stressMagnitude1 - 2.9645), epsilon);

    auto stress2 = wind.calculate(Wind::SpeedVector {
        .u10 = 22,
        .v10 = 0
    });

    BOOST_CHECK_EQUAL(stress2.qx, 1.48225);
    BOOST_CHECK_EQUAL(stress2.qy, 0);
}
#include <format>

#include <doctest/doctest.h>

#include <calc/turbulence.h>

using namespace std;

using namespace Calc;

TEST_CASE("Test default viscosity params") {
    Turbulence::DefaultViscosityParams params {
        .turbulence = Turbulence::Params {
            .k0 = 0.4,
            .k = 2e-5,
            .sigma = 1.2e-4,
            .num = 0.0001
        },
        .phys = Data::PhysicalParams {
            .rho = 1025,
            .f = 1.2e-4
        }
    };

    params.turbulence.k0 = 0;

    CHECK_THROWS_WITH_AS(
        Turbulence::DefaultViscosity viscosity(params),
        format("k0 should be > 0, but it is {}", params.turbulence.k0).c_str(),
        runtime_error
    );    

    params.turbulence.k0 = 0.4;
    params.turbulence.k = 0;

    CHECK_THROWS_WITH_AS(
        Turbulence::DefaultViscosity viscosity(params),
        format("k should be > 0, but it is {}", params.turbulence.k).c_str(),
        runtime_error
    );

    params.turbulence.k = 2e-5;
    params.turbulence.sigma = 0;

    CHECK_THROWS_WITH_AS(
        Turbulence::DefaultViscosity viscosity(params),
        format("sigma should be > 0, but it is {}", params.turbulence.sigma).c_str(),
        runtime_error
    );

    params.turbulence.sigma = 1.2e-4;
    params.turbulence.num = 0;

    CHECK_THROWS_WITH_AS(
        Turbulence::DefaultViscosity viscosity(params),
        format("nu should be > 0, but it is {}", params.turbulence.num).c_str(),
        runtime_error
    );

    params.turbulence.num = 0.0001;
    params.phys.rho = 0;

    CHECK_THROWS_WITH_AS(
        Turbulence::DefaultViscosity viscosity(params),
        format("rho should be > 0, but it is {}", params.phys.rho).c_str(),
        runtime_error
    );

    params.phys.rho = 1025;
    params.phys.f = 0;

    CHECK_THROWS_WITH_AS(
        Turbulence::DefaultViscosity viscosity(params),
        format("f should be > 0, but it is {}", params.phys.f).c_str(),
        runtime_error
    );
}

TEST_CASE("Test default viscosity calc") {
    double epsilon = 1e-15;

    double num = 0.0001;

    Turbulence::DefaultViscosityParams params {
        .turbulence = Turbulence::Params {
            .k0 = 0.4,
            .k = 2e-5,
            .sigma = 1.2e-4,
            .num = num
        },
        .phys = Data::PhysicalParams {
            .rho = 1025,
            .f = 1.2e-4
        }
    };

    Turbulence::DefaultViscosity viscosity(params);

    Turbulence::WindCurrentParams windCurrent {
        .wind = Wind::SpeedVector {
            .u10 = 22,
            .v10 = 0
        },
        .stress = Wind::StressVector {
            .qx = 1.48225,
            .qy = 0
        },
        .current = Data::DepthAveragedVelocity {
            .ua = 0.5,
            .va = 0
        }               
    };

    auto fricInf = viscosity.calculate(windCurrent);

    CHECK_LT(abs(fricInf.ut - 0.038027589471009204), epsilon);    
    CHECK_LT(abs(fricInf.ht - 126.75863157003069), epsilon);

    auto friction = Turbulence::FrictionParams {
        .fricInf = fricInf,
        .windCurrent = windCurrent
    };

    auto viscosityState = viscosity.calculate(friction);

    CHECK_LT(abs(viscosityState.nut - 0.041666666666666671 - num), epsilon);
    CHECK_LT(abs(viscosityState.nus - 0.41430442406520351 - num), epsilon);

    Turbulence::ViscosityParams viscosityParams {
        .viscosity = viscosityState,
        .friction = friction
    };

    CHECK_EQ(
        viscosity.calculate(0, viscosityParams),
        viscosityState.nus
    );

    CHECK_EQ(
        viscosity.calculate(fricInf.ht + epsilon, viscosityParams),
        viscosityState.nut
    );
}
#include <calc/turbulence.h>

using namespace Calc;

BOOST_AUTO_TEST_CASE(test_default_viscosity_params) {
    Turbulence::DefaultViscosityParams params {
        .turbulence = Turbulence::TurbulenceParams {
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

    BOOST_CHECK_EXCEPTION(
        Turbulence::DefaultViscosity viscosity(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("k0 should be > 0, but it is {}", params.turbulence.k0);
        }
    );    

    params.turbulence.k0 = 0.4;
    params.turbulence.k = 0;

    BOOST_CHECK_EXCEPTION(
        Turbulence::DefaultViscosity viscosity(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("k should be > 0, but it is {}", params.turbulence.k);
        }
    );

    params.turbulence.k = 2e-5;
    params.turbulence.sigma = 0;

    BOOST_CHECK_EXCEPTION(
        Turbulence::DefaultViscosity viscosity(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("sigma should be > 0, but it is {}", params.turbulence.sigma);
        }
    );

    params.turbulence.sigma = 1.2e-4;
    params.turbulence.num = 0;

    BOOST_CHECK_EXCEPTION(
        Turbulence::DefaultViscosity viscosity(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("nu should be > 0, but it is {}", params.turbulence.num);
        }
    );

    params.turbulence.num = 0.0001;
    params.phys.rho = 0;

    BOOST_CHECK_EXCEPTION(
        Turbulence::DefaultViscosity viscosity(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("rho should be > 0, but it is {}", params.phys.rho);
        }
    );

    params.phys.rho = 1025;
    params.phys.f = 0;

    BOOST_CHECK_EXCEPTION(
        Turbulence::DefaultViscosity viscosity(params),
        runtime_error,
        [params](const runtime_error& e) {
            return string(e.what()) == format("f should be > 0, but it is {}", params.phys.f);
        }
    );
}

BOOST_AUTO_TEST_CASE(test_default_viscosity_calc) {
    double epsilon = 1e-15;

    double num = 0.0001;

    Turbulence::DefaultViscosityParams params {
        .turbulence = Turbulence::TurbulenceParams {
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

    BOOST_CHECK_LT(abs(fricInf.ut - 0.038027589471009204), epsilon);    
    BOOST_CHECK_LT(abs(fricInf.ht - 126.75863157003069), epsilon);

    auto friction = Turbulence::FrictionParams {
        .fricInf = fricInf,
        .windCurrent = windCurrent
    };

    auto viscosityState = viscosity.calculate(friction);

    BOOST_CHECK_LT(abs(viscosityState.nut - 0.041666666666666671 - num), epsilon);
    BOOST_CHECK_LT(abs(viscosityState.nus - 0.41430442406520351 - num), epsilon);

    Turbulence::ViscosityParams viscosityParams {
        .viscosity = viscosityState,
        .friction = friction
    };

    BOOST_CHECK_EQUAL(
        viscosity.calculate(0, viscosityParams),
        viscosityState.nus
    );

    BOOST_CHECK_EQUAL(
        viscosity.calculate(fricInf.ht + epsilon, viscosityParams),
        viscosityState.nut
    );
}
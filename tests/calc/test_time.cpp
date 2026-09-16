#include <calc/time.h>

using namespace Calc;

BOOST_AUTO_TEST_CASE(test_default_adjust_time_step_params) {
	Time::DefaultAdjustTimeStepParams params{
		.b = -1,
		.tMax = -1,
		.dtMax = -1
	};

	BOOST_CHECK_EXCEPTION(
		Time::DefaultAdjustTimeStep ts(params),
		runtime_error,
		[params](const runtime_error& e) {
			return string(e.what()) == format("b should be > 0, but it is {}", params.b);
		}
	);

	params.b = 1;

	BOOST_CHECK_EXCEPTION(
		Time::DefaultAdjustTimeStep ts(params),
		runtime_error,
		[params](const runtime_error& e) {
			return string(e.what()) == format("tMax should be > 0, but it is {}", params.tMax);
		}
	);

	params.tMax = 1;

	BOOST_CHECK_EXCEPTION(
		Time::DefaultAdjustTimeStep ts(params),
		runtime_error,
		[params](const runtime_error& e) {
			return string(e.what()) == format("dtMax should be > 0, but it is {}", params.dtMax);
		}
	);
}

BOOST_AUTO_TEST_CASE(test_default_adjust_time_step_calc) {
	Time::DefaultAdjustTimeStepParams params {
		.b = 1.1,
		.tMax = 0.3,
		.dtMax = 0.03
	};

	Time::DefaultAdjustTimeStep ts(params);

	Time::State state1 {
		.t = 0,
		.dt = 0.04,
		.mult = false
	};

	BOOST_CHECK_EQUAL(ts.calculate(state1), params.dtMax);

	Time::State state2 {
		.t = 0.28,
		.dt = 0.029,
		.mult = false
	};

	BOOST_CHECK_EQUAL(ts.calculate(state2), params.tMax - state2.t);

	Time::State state3 {
		.t = 0,
		.dt = 0.01,
		.mult = true
	};

	BOOST_CHECK_EQUAL(ts.calculate(state3), params.b * state3.dt);

	Time::State state4 {
		.t = 0,
		.dt = 0.01,
		.mult = false
	};

	BOOST_CHECK_EQUAL(ts.calculate(state4), state4.dt);
}
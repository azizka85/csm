#include <format>

#include <doctest/doctest.h>

#include <calc/time.h>

using namespace std;

using namespace Calc;

TEST_CASE("Test default adjust time step params") {
	Time::DefaultAdjustTimeStepParams params{
		.b = -1,
		.tMax = -1
	};

	CHECK_THROWS_WITH_AS(
		Time::DefaultAdjustTimeStep ts(params),
		format("b should be > 0, but it is {}", params.b).c_str(),
		runtime_error
	);

	params.b = 1;

	CHECK_THROWS_WITH_AS(
		Time::DefaultAdjustTimeStep ts(params),
		format("tMax should be > 0, but it is {}", params.tMax).c_str(),
		runtime_error
	);
}

TEST_CASE("Test default adjust time step calc") {
	Time::DefaultAdjustTimeStepParams params {
		.b = 1.1,
		.tMax = 0.3		
	};

	Time::DefaultAdjustTimeStep ts(params);

	Time::State state1 {
		.t = 0,
		.dt = 0.04,
		.dtMax = 0.03,
		.mult = false
	};

	CHECK_EQ(ts.calculate(state1), state1.dtMax);

	Time::State state2 {
		.t = 0.28,
		.dt = 0.029,
		.dtMax = 0.03,
		.mult = false
	};

	CHECK_EQ(ts.calculate(state2), params.tMax - state2.t);

	Time::State state3 {
		.t = 0,
		.dt = 0.01,
		.dtMax = 0.03,
		.mult = true
	};

	CHECK_EQ(ts.calculate(state3), params.b * state3.dt);

	Time::State state4 {
		.t = 0,
		.dt = 0.01,
		.dtMax = 0.03,
		.mult = false
	};

	CHECK_EQ(ts.calculate(state4), state4.dt);
}
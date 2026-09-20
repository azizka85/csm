#include <format>

#include <utils/checks/time.h>

#include "time.h"

using namespace std;

using namespace Utils;

Calc::Time::DefaultAdjustTimeStep::DefaultAdjustTimeStep(DefaultAdjustTimeStepParams params) {
	Checks::Time::checkB(params.b);
	Checks::Time::checkTMax(params.tMax);

	this->params = move(params);
}

double Calc::Time::DefaultAdjustTimeStep::calculate(Calc::Time::State state) {
	auto b = params.b;
	auto tMax = params.tMax;	

	auto t = state.t;
	auto dt = state.dt;
	auto dtMax = state.dtMax;
	auto mult = state.mult;

	if (dt >= dtMax) {
		return dtMax;
	}

	if (t < tMax && t + dt >= tMax) {
		return tMax - t;
	}

	if (mult) {
		return b * dt;
	}

	return dt;
}

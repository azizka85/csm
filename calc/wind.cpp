#include <cmath>

#include <format>

#include <utils/checks/phys.h>

#include "wind.h"

using namespace Utils::Checks;
using namespace Calc;

Wind::DefaultStress::DefaultStress(DefaultStressParams params) {
	Phys::checkRho(params.rhoAir);
	Phys::checkCD(params.cd);

	this->params = move(params);
}

Wind::StressVector Wind::DefaultStress::calculate(SpeedVector windSpeed) {
	double rhoAir = params.rhoAir;
	double cd = params.cd;

	double u10 = windSpeed.u10;
	double v10 = windSpeed.v10;

	double magnitude = sqrt(u10 * u10 + v10 * v10);

	return StressVector {
		.qx = cd * rhoAir * u10 * magnitude,
		.qy = cd * rhoAir * v10 * magnitude
	};
}

string Wind::DefaultStress::dirName(string attr) {
	return format("{}, rho_air={}, cd={}", attr, params.rhoAir, params.cd);
}

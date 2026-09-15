#include <format>

#include <utils/checks/bathymetry.h>

#include "grid.h"

using namespace std;

using namespace Calc::Grid;
using namespace Utils::Checks;

TriplePointVerticalStepSize::TriplePointVerticalStepSize(TriplePointVerticalStepSizeParams params) {
	Bathymetry::checkZM(params.zm);
	Bathymetry::checkDZMinMax(params.dzMin, params.dzMax);

	auto zm = params.zm;
	auto dzMax = params.dzMax;

	this->k =
		dzMax / (zm * zm * zm * zm / 4 - (zm + 1) * zm * zm * zm / 3 + zm * zm * zm / 2);

	this->params = move(params);
}

string TriplePointVerticalStepSize::dirName(string attr) {
	return format("{}, dz={}-{}, z={}", attr, params.dzMin, params.dzMax, params.zm);
}

double TriplePointVerticalStepSize::calculate(double z) {
	auto zm = params.zm;
	auto dzMin = params.dzMin;

	return k * (z * z * z * z / 4 - (zm + 1) * z * z * z / 3 + zm * z * z / 2) + dzMin;
}

#include <format>

#include <utils/checks/phys.h>

#include "turbulence.h"

using namespace Utils::Checks;

using namespace Calc;

Turbulence::DefaultViscosity::DefaultViscosity(DefaultViscosityParams params) {
	auto& turbulence = params.turbulence;
	auto& phys = params.phys;

	Phys::checkK0(turbulence.k0);
	Phys::checkK(turbulence.k);
	Phys::checkSigma(turbulence.sigma);
	Phys::checkNU(turbulence.num);

	Phys::checkRho(phys.rho);
	Phys::checkF(phys.f);

	this->params = move(params);
}

Turbulence::FrictionalInfluenceParams Turbulence::DefaultViscosity::calculate(WindCurrentParams windCurrent) {
	auto& turbulence = params.turbulence;
	auto& phys = params.phys;

	auto& stress = windCurrent.stress;

	double k0 = turbulence.k0;
	double rho = phys.rho;
	double f = phys.f;

	double qx = stress.qx;
	double qy = stress.qy;

	auto q = sqrt(qx * qx + qy * qy);

	auto ut = sqrt(q / rho);
	auto ht = k0 * ut / f;

	return FrictionalInfluenceParams {
		.ut = ut,
		.ht = ht
	};
}

Turbulence::ViscosityState Turbulence::DefaultViscosity::calculate(FrictionParams friction) {
	auto& turbulence = params.turbulence;

	auto& wind = friction.windCurrent.wind;
	auto& current = friction.windCurrent.current;

	auto u10 = wind.u10;
	auto v10 = wind.v10;

	double ua = current.ua;
	double va = current.va;

	double k = turbulence.k;
	double sigma = turbulence.sigma;
	double num = turbulence.num;

	auto nut = k * (ua * ua + va * va) / sigma + num;

	auto wind_mag_sq = u10 * u10 + v10 * v10;

	auto nus = 0.1825e-3 * pow(wind_mag_sq, 1.25) + num;

	return ViscosityState {
		.nut = nut,
		.nus = nus
	};
}

string Turbulence::DefaultViscosity::dirName(string attr) {
	auto& turbulence = params.turbulence;

	return format(
		"{}, k0={}, k={}, sigma={}, nu0={}",
		attr, 
		turbulence.k0, turbulence.k, turbulence.sigma, turbulence.num
	);
}

double Turbulence::DefaultViscosity::calculate(double z, ViscosityParams viscosity) {
	double ht = viscosity.friction.fricInf.ht;

	double nus = viscosity.viscosity.nus;
	double nut = viscosity.viscosity.nut;

	if (z < ht) {
		return nus - (nus - nut) * z / ht;
	}
	
	return nut;
}

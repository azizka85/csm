#ifndef CALC_TURBULENCE_H
#define CALC_TURBULENCE_H

#include "data.h"
#include "wind.h"

using namespace Calc;

namespace Calc::Turbulence {
	struct FrictionalInfluenceParams {
		double ut;
		double ht;
	};

	struct ViscosityState {
		double nut;
		double nus;
	};

	struct WindCurrentParams {
		Wind::SpeedVector wind;
		Wind::StressVector stress;
		Data::DepthAveragedVelocity current;
	};

	struct FrictionParams {
		FrictionalInfluenceParams fricInf;
		WindCurrentParams windCurrent;
	};

	struct ViscosityParams {
		ViscosityState viscosity;
		FrictionParams friction;
	};

	struct TurbulenceParams {
		double k0;
		double k;
		double sigma;
		double num;
	};

	class IFrictionalInfluence {
		public:
			virtual FrictionalInfluenceParams calculate(WindCurrentParams params) = 0;
	};

	class ITidalSurfaceViscosity {
		public:
			virtual ViscosityState calculate(FrictionParams params) = 0;
	};

	class IViscosity {
		public:
			virtual string dirName(string attr) = 0;
			virtual double calculate(double z, ViscosityParams params) = 0;
	};

	struct DefaultViscosityParams {
		TurbulenceParams turbulence;
		Data::PhysicalParams phys;
	};

	class DefaultViscosity : IFrictionalInfluence, ITidalSurfaceViscosity, IViscosity {
		private:
			DefaultViscosityParams params;

		public:
			DefaultViscosity(DefaultViscosityParams params);

			FrictionalInfluenceParams calculate(WindCurrentParams windCurrent) override;

			ViscosityState calculate(FrictionParams friction) override;

			string dirName(string attr) override;
			double calculate(double z, ViscosityParams viscosity) override;
	};
}

#endif
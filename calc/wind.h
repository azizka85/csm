#ifndef CALC_WIND_H
#define CALC_WIND_H

#include <string>

using namespace std;

namespace Calc::Wind {
	struct SpeedVector {
		double u10;
		double v10;
	};

	struct StressVector {
		double qx;
		double qy;
	};

	class IStress {
		public:
			virtual StressVector calculate(SpeedVector windSpeed) = 0;
			virtual string dirName(string attr) = 0;
	};

	struct DefaultStressParams {
		double rhoAir;
		double cd;
	};

	class DefaultStress: IStress {
		private:
			DefaultStressParams params;

		public:	
			DefaultStress(DefaultStressParams params);

			StressVector calculate(SpeedVector windSpeed) override;
			string dirName(string attr) override;

	};
}

#endif
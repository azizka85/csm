#ifndef CALC_GRID_H
#define CALC_GRID_H

#include <string>

using namespace std;

namespace Calc::Grid {
	class INonUniformVerticalStepSize {
		public:
			virtual string dirName(string attr) = 0;
			virtual double calculate(double z) = 0;
	};

	struct TriplePointVerticalStepSizeParams {
		double zm;
		double dzMin;
		double dzMax;
	};

	class TriplePointVerticalStepSize : public INonUniformVerticalStepSize {
		private:
			TriplePointVerticalStepSizeParams params;

			double k;

		public:
			TriplePointVerticalStepSize(TriplePointVerticalStepSizeParams params);

			string dirName(string attr) override;
			double calculate(double z) override;

	};
}

#endif
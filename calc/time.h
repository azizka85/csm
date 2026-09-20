#ifndef CALC_TIME_H
#define CALC_TIME_H

namespace Calc::Time {
	struct State {
		double t;
		double dt;
		double dtMax;
		bool mult;
	};

	class IAdjustTimeStep {
		public:
			virtual double calculate(State state) = 0;
	};

	struct DefaultAdjustTimeStepParams {
		double b;
		double tMax;				
	};

	class DefaultAdjustTimeStep : public IAdjustTimeStep {
		private:
			DefaultAdjustTimeStepParams params;

		public:		
			DefaultAdjustTimeStep(DefaultAdjustTimeStepParams params);

			double calculate(State state) override;

	};
}

#endif
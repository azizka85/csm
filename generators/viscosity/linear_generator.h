#ifndef GENERATORS_VISCOSITY_LINEAR_GENERATOR_H
#define GENERATORS_VISCOSITY_LINEAR_GENERATOR_H

#include <calc/turbulence.h>

#include "../viscosity.h"



namespace Generators::Viscosity {
	class LinearGenerator : IGenerator {
		private:
			double ht;
			Turbulence::ViscosityState viscosity;

		public:
			LinearGenerator(double ht, Turbulence::ViscosityState viscosity);

			path addDirectory(path outDir) override;
			vector<double> generate(Data::VolumeState volume, Data::WindCurrent windCurrent) override;
	};
}

#endif 
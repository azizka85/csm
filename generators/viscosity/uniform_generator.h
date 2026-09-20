#ifndef GENERATORS_VISCOSITY_UNIFORM_GENERATOR_H
#define GENERATORS_VISCOSITY_UNIFORM_GENERATOR_H

#include "../viscosity.h"

namespace Generators::Viscosity {
	class UniformGenerator: public IGenerator {
		private:
			double nu;

		public:
			UniformGenerator(double nu);

			path addDirectory(path outDir) override;
			vector<double> generate(Data::VolumeState volume, Data::WindCurrent windCurrent) override;
	};
}

#endif 
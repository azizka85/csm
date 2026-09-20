#ifndef GENERATORS_DZ_UNIFORM_GENERATOR_H
#define GENERATORS_DZ_UNIFORM_GENERATOR_H

#include "../dz.h"

namespace Generators::DZ {
	class UniformGenerator: public IGenerator {
		private:
			double dz;

		public:
			UniformGenerator(double dz);

			path addDirectory(path outDir) override;
			vector<double> generate() override;
	};
}

#endif 
#ifndef GENERATORS_WIND_UNIFORM_GENERATOR_H
#define GENERATORS_WIND_UNIFORM_GENERATOR_H

#include <calc/wind.h>

#include "../wind.h"

using namespace std::filesystem;

namespace Generators::Wind {
	struct UniformGeneratorParams {
		Calc::Wind::SpeedVector speed;
		Calc::Wind::StressVector stress;
	};

	class UniformGenerator: IGenerator {
		private:
			UniformGeneratorParams params;

		public:
			UniformGenerator(UniformGeneratorParams params);

			virtual path addDirectory(path outDir) = 0;
			virtual vector<Data::Wind> generate(size_t nx, size_t ny) = 0;
	};
}

#endif 
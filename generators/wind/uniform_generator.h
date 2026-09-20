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

	class UniformGenerator: public IGenerator {
		private:
			UniformGeneratorParams params;

		public:
			UniformGenerator(UniformGeneratorParams params);
			path addDirectory(path outDir) override;
			vector<Data::Wind> generate(size_t nx, size_t ny) override;
	};
}

#endif 
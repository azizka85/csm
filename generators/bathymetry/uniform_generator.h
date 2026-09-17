#ifndef GENERATORS_BATHYMETRY_UNIFORM_GENERATOR_H
#define GENERATORS_BATHYMETRY_UNIFORM_GENERATOR_H

#include <filesystem>

#include "../bathymetry.h"

using namespace std::filesystem;

namespace Generators::Bathymetry {
	class UniformGenerator : IGenerator {
		private:
			double h;

		public:
			UniformGenerator(double h);			

			path addDirectory(path outDir) override;
			vector<double> generate(double nx, double ny) override;
	};
}

#endif 
#ifndef GENERATORS_BATHYMETRY_COSINE_GENERATOR_H
#define GENERATORS_BATHYMETRY_COSINE_GENERATOR_H

#include <filesystem>

#include "../bathymetry.h"

using namespace std::filesystem;

namespace Generators::Bathymetry {
	class CosineGenerator : public IGenerator {
		private:
			double hm;

		public:
			CosineGenerator(double hm);

			path addDirectory(path outDir) override;
			vector<double> generate(double nx, double ny) override;
	};
}

#endif 
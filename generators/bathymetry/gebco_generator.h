#ifndef GENERATORS_BATHYMETRY_GEBCO_GENERATOR_H
#define GENERATORS_BATHYMETRY_GEBCO_GENERATOR_H

#include <filesystem>

#include <calc/data.h>

#include "../bathymetry.h"

using namespace std::filesystem;

using namespace Calc;

namespace Generators::Bathymetry {
	struct GEBCOGeneratorParams {
		Data::GeoBound bound;
		double refDepth;
		double minDepth;
		path filePath;
	};

	class GEBCOGenerator : public IGenerator {
		private:
			GEBCOGeneratorParams params;

		public:
			GEBCOGenerator(GEBCOGeneratorParams params);

			path addDirectory(path outDir) override;
			vector<double> generate(double nx, double ny) override;
	};
}

#endif
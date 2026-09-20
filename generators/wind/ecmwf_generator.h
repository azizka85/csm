#ifndef GENERATORS_WIND_ECMWF_GENERATOR_H
#define GENERATORS_WIND_ECMWF_GENERATOR_H

#include <calc/wind.h>

#include "../wind.h"

using namespace std::filesystem;

using namespace Calc::Wind;

namespace Generators::Wind {
	struct ECMWFGeneratorParams {
		Data::GeoBound bound;
		shared_ptr<IStress> windStress;
		path filePath;
		string dirAttr;
	};

	class ECMWFGenerator : public IGenerator {
		private:
			ECMWFGeneratorParams params;

		public:
			ECMWFGenerator(ECMWFGeneratorParams params);

			path addDirectory(path outDir) override;
			vector<Data::Wind> generate(size_t nx, size_t ny) override;
	};
}

#endif 
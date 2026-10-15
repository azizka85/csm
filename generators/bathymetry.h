#ifndef GENERATORS_BATHYMETRY_H
#define GENERATORS_BATHYMETRY_H

#include <vector>

#include <filesystem>

#include <calc/data.h>

using namespace std;
using namespace std::filesystem;

using namespace Calc;

namespace Generators::Bathymetry {
	class IGenerator {
		public:
			virtual path addDirectory(path outDir) = 0;
			
			virtual Data::Bathymetry generate(
				size_t nx, size_t ny,
				double dx, double dy
			) = 0;
	};
}

#endif 
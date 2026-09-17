#ifndef GENERATORS_BATHYMETRY_H
#define GENERATORS_BATHYMETRY_H

#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace Generators::Bathymetry {
	class IGenerator {
		public:
			virtual path addDirectory(path outDir) = 0;
			virtual vector<double> generate(double nx, double ny) = 0;
	};
}

#endif 
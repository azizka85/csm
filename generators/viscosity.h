#ifndef GENERATORS_VISCOSITY_H
#define GENERATORS_VISCOSITY_H

#include <vector>

#include <filesystem>

#include <calc/data.h>

using namespace std;
using namespace std::filesystem;

using namespace Calc;

namespace Generators::Viscosity {
	class IGenerator {
		public:
			virtual path addDirectory(path outDir) = 0;
			virtual vector<double> generate(Data::VolumeState volume, Data::WindCurrent windCurrent) = 0;
	};
}

#endif 
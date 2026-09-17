#ifndef GENERATORS_WIND_H
#define GENERATORS_WIND_H

#include <filesystem>

#include <calc/data.h>

using namespace std::filesystem;

using namespace Calc;

namespace Generators::Wind {
	class IGenerator {
		public:
			virtual path addDirectory(path outDir) = 0;
			virtual vector<Data::Wind> generate(size_t nx, size_t ny) = 0;
	};
}

#endif 
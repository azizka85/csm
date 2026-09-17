#ifndef GENERATORS_DZ_H
#define GENERATORS_DZ_H

#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace Generators::DZ {
	class IGenerator {
		public:
			virtual path addDirectory(path outDir) = 0;
			virtual vector<double> generate() = 0;
	};
}

#endif 
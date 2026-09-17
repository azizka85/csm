#ifndef GENERATORS_AREA_H
#define GENERATORS_AREA_H

#include <filesystem>

using namespace std::filesystem;

namespace Generators::Area {
	struct Geometry {
		double l;
		double w;
	};

	class IGenerator {
		public:
			virtual path addDirectory(path outDir) = 0;
			virtual Geometry generate() = 0;
	};
}

#endif 
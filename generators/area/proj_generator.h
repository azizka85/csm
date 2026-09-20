#ifndef GENERATORS_AREA_PROJ_GENERATOR_H
#define GENERATORS_AREA_PROJ_GENERATOR_H

#include <calc/data.h>

#include "../area.h"

using namespace Calc;

namespace Generators::Area {	
	class ProjGenerator : public IGenerator {
		private:
			Data::GeoBound bound;

		public:
			ProjGenerator(Data::GeoBound bound);

			path addDirectory(path outDir) override;
			Geometry generate() override;
	};
}

#endif 
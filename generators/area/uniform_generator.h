#ifndef GENERATORS_AREA_UNIFORM_GENERATOR_H
#define GENERATORS_AREA_UNIFORM_GENERATOR_H

#include "../area.h"

namespace Generators::Area {
	class UniformGenerator : public IGenerator {
		private:
			Geometry geom;

		public:
			UniformGenerator(Geometry geom);

			path addDirectory(path outDir) override;
			Geometry generate() override;

	};
}

#endif 
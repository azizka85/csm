#ifndef GENERATORS_DZ_TRIPLE_POINT_GENERATOR_H
#define GENERATORS_DZ_TRIPLE_POINT_GENERATOR_H

#include <calc/grid.h>

#include "../dz.h"

using namespace Calc;

namespace Generators::DZ {
	struct TriplePointGeneratorParams {
		shared_ptr<Grid::INonUniformVerticalStepSize> vertStep;

		string dirAttr;
	};

	class TriplePointGenerator : public IGenerator {
		private:
			TriplePointGeneratorParams params;

		public:
			TriplePointGenerator(TriplePointGeneratorParams params);

			path addDirectory(path outDir) override;
			vector<double> generate() override;
	};
}

#endif 
#include "uniform_generator.h"

using namespace Generators::Wind;

UniformGenerator::UniformGenerator(UniformGeneratorParams params) {
	this->params = move(params);
}

path UniformGenerator::addDirectory(path outDir) {
	auto qx = params.stress.qx;
	auto qy = params.stress.qy;

	return outDir / path(
		format("uq, qx={}, qy={}", qx, qy)
	);
}

vector<Data::Wind> UniformGenerator::generate(size_t nx, size_t ny) {
	return {
		Data::Wind {
			.time = 0,
			.speed = Data::WindSpeed {
				.u10 = vector<double>(nx * ny, params.speed.u10),
				.v10 = vector<double>(nx * ny, params.speed.v10)
			},
			.stress = Data::WindStress {
				.qx = vector<double>(nx * ny, params.stress.qx),
				.qy = vector<double>(nx * ny, params.stress.qy)
			}
		}
	};
}
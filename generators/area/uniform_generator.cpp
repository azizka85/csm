#include <format>

#include <utils/checks/area.h>

#include "uniform_generator.h"
#include "uniform_generator.h"
#include "uniform_generator.h"

using namespace std;

using namespace Utils;

using namespace Generators;

Area::UniformGenerator::UniformGenerator(Geometry geom) {
	Checks::Area::checkL(geom.l);
	Checks::Area::checkW(geom.w);

	this->geom = move(geom);
}

path Area::UniformGenerator::addDirectory(path outDir) {
	return outDir / path(
		format("l={}, w={}", geom.l, geom.w)
	);
}

Area::Geometry Area::UniformGenerator::generate() {
	return geom;
}

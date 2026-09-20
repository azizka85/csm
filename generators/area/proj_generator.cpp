#include <utils/checks/area.h>

#include "../utils.h"

#include "proj_generator.h"

Generators::Area::ProjGenerator::ProjGenerator(Data::GeoBound bound) {
	::Utils::Checks::Area::checkLatMinMax(bound.latMin, bound.latMax);
	::Utils::Checks::Area::checkLonMinMax(bound.lonMin, bound.lonMax);

	this->bound = std::move(bound);
}

path Generators::Area::ProjGenerator::addDirectory(path outDir) {
	return outDir / path(
		format("lat={}-{}, lon={}-{}", bound.latMin, bound.latMax, bound.lonMin, bound.lonMax)
	);
}

Generators::Area::Geometry Generators::Area::ProjGenerator::generate() {
	auto cartBound = Utils::toCartesian(bound);

	return Geometry{
		.l = cartBound.maxX - cartBound.minX,
		.w = cartBound.maxY - cartBound.minY
	};
}

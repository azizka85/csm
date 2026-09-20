#include <doctest/doctest.h>

#include <generators/utils.h>

using namespace Calc::Data;

using namespace Generators::Utils;

TEST_CASE("Test geo to cartesian") {
	double epsilon = 1e-8;

	auto bound = toCartesian(GeoBound{
		.latMin = 36,
		.lonMin = 46,
		
		.latMax = 47,
		.lonMax = 56
	});

	CHECK_LT(abs(bound.minX - 49183.89173453045), epsilon);
	CHECK_LT(abs(bound.maxX - 950816.108265473), epsilon);

	CHECK_LT(abs(bound.minY - 3995527.639948009), epsilon);
	CHECK_LT(abs(bound.maxY - 5217308.379028135), epsilon);
}

TEST_CASE("Test cartesian to geo") {
	double epsilon = 1e-8;

	auto bound = toGeo(CartesianBound{
		.minX = 49183.89173453045,
		.minY = 3995527.639948009,

		.maxX = 950816.108265473,
		.maxY = 5217308.379028135
	});

	CHECK_LT(abs(bound.latMin - 36), epsilon);
	CHECK_LT(abs(bound.latMax - 46.955704812277865), epsilon);

	CHECK_LT(abs(bound.lonMin - 45.074376577489787), epsilon);
	CHECK_LT(abs(bound.lonMax - 56.925623422510263), epsilon);
}
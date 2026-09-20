#include <format>

#include <doctest/doctest.h>

#include <generators/area/proj_generator.h>

using namespace Generators::Area;

TEST_CASE("Test area proj generator params") {
	double epsilon = 1e-8;

	Data::GeoBound bound{
		.latMin = 36,
		.lonMin = 46,

		.latMax = 27,
		.lonMax = 56
	};

	CHECK_THROWS_WITH_AS(
		ProjGenerator generator(bound),
		format("latMin should be < latMax, but they are {} and {}", bound.latMin, bound.latMax).c_str(),
		runtime_error
	);

	bound = Data::GeoBound{
		.latMin = 36,
		.lonMin = 46,

		.latMax = 47,
		.lonMax = 46
	};

	CHECK_THROWS_WITH_AS(
		ProjGenerator generator(bound),
		format("lonMin should be < lonMax, but they are {} and {}", bound.lonMin, bound.lonMax).c_str(),
		runtime_error
	);
}

TEST_CASE("Test area proj generator generate") {
	double epsilon = 1e-8;

	Data::GeoBound bound{
		.latMin = 36,
		.lonMin = 46,

		.latMax = 47,
		.lonMax = 56
	};

	auto geom = ProjGenerator(bound).generate();

	CHECK_LT(abs(geom.l - 901632.21653094255), epsilon);
	CHECK_LT(abs(geom.w - 1221780.739080126), epsilon);
}
#include <format>

#include <doctest/doctest.h>

#include <generators/area/uniform_generator.h>

using namespace std;

using namespace Generators::Area;

TEST_CASE("Test area uniform generator params") {
	Geometry geom { .l = -1, .w = 5 };

	CHECK_THROWS_WITH_AS(
		UniformGenerator area(geom),
		format("l should be > 0, but it is {}", geom.l).c_str(),
		runtime_error
	);

	geom = Geometry { .l = 5, .w = -1 };

	CHECK_THROWS_WITH_AS(
		UniformGenerator area(geom),
		format("w should be > 0, but it is {}", geom.w).c_str(),
		runtime_error
	);
}

TEST_CASE("Test area uniform generator generate") {
	Geometry geomIn { .l = 7.2, .w = 6.7 };
	
	UniformGenerator generator(geomIn);

	auto geomOut = generator.generate();

	CHECK_EQ(geomOut.l, geomIn.l);
	CHECK_EQ(geomOut.w, geomIn.w);
}
#include <doctest/doctest.h>

#include <generators/bathymetry/uniform_generator.h>

using namespace Generators::Bathymetry;

TEST_CASE("Test bathymetry uniform generator params") {
	double h = -35;

	CHECK_THROWS_WITH_AS(
		UniformGenerator generator(h),
		format("h should be > 0, but it is {}", h).c_str(),
		runtime_error
	);	
}

TEST_CASE("Test bathymetry uniform generator generate") {
	size_t nx = 1;
	size_t ny = 2;
	
	double h = 260;
	
	UniformGenerator generator(h);

	auto res = generator.generate(nx, ny);

	CHECK_EQ(res.size(), nx * ny);

	for (size_t i = 0; i < nx * ny; i++) {
		CHECK_EQ(res[i], h);
	}
}
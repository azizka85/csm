#include <doctest/doctest.h>

#include <generators/dz/uniform_generator.h>

using namespace Generators::DZ;

TEST_CASE("Test dz uniform generator params") {
	double dz = 0;

	CHECK_THROWS_WITH_AS(
		UniformGenerator generator(dz),
		format("dz should be > 0, but it is {}", dz).c_str(),
		runtime_error
	);
}

TEST_CASE("Test dz uniform generator generate") {
	double dz = 0.002;

	size_t n = ceil(1 / dz);

	UniformGenerator generator(dz);

	auto res = generator.generate();

	CHECK_EQ(res.size(), n);

	for (size_t i = 0; i < n; i++) {
		CHECK_EQ(res[i], dz);
	}
}
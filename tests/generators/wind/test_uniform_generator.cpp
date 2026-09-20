#include <doctest/doctest.h>

#include <generators/wind/uniform_generator.h>

using namespace Calc::Wind;
using namespace Generators::Wind;

TEST_CASE("Test wind uniform generator generate") {
	SpeedVector speed{ .u10 = 22, .v10 = 22 };
	StressVector stress{ .qx = 1.5, .qy = 1.5 };

	UniformGenerator generator(
		UniformGeneratorParams{
			.speed = speed,
			.stress = stress
		}
	);

	size_t nx = 11;
	size_t ny = 11;

	auto data = generator.generate(nx, ny);

	CHECK_EQ(data.size(), 1);

	CHECK_EQ(data[0].time, 0);

	CHECK_EQ(data[0].speed.u10.size(), nx * ny);
	CHECK_EQ(data[0].speed.v10.size(), nx * ny);

	CHECK_EQ(data[0].stress.qx.size(), nx * ny);
	CHECK_EQ(data[0].stress.qy.size(), nx * ny);

	for (size_t i = 0; i < nx*ny; i++) {
		CHECK_EQ(data[0].speed.u10[i], speed.u10);
		CHECK_EQ(data[0].speed.v10[i], speed.v10);

		CHECK_EQ(data[0].stress.qx[i], stress.qx);
		CHECK_EQ(data[0].stress.qy[i], stress.qy);
	}
}
#include <doctest/doctest.h>

#include <generators/viscosity/uniform_generator.h>

using namespace Calc::Data;
using namespace Generators::Viscosity;

TEST_CASE("Test viscosity uniform generator params") {
	double nu = 0;

	CHECK_THROWS_WITH_AS(
		UniformGenerator generator(nu),
		format("nu should be > 0, but it is {}", nu).c_str(),
		runtime_error
	);
}

TEST_CASE("Test viscosity uniform generator generate") {
	size_t nx = 10;
	size_t ny = 11;
	size_t nz = 12;

	double nu = 0.4;

	UniformGenerator generator(nu);

	vector<double> v;

	auto res = generator.generate(
		VolumeState{
			.nx = nx,
			.ny = ny,
			.nz = nz,
			.dx = 0,
			.dy = 0,
			.dz = v,
			.h = v
		},
		WindCurrent{
			.speed = WindSpeed {
				.u10 = {},
				.v10 = {}
			},
			.stress = WindStress {
				.qx = {},
				.qy = {}
			},
			.current = Current {
				.ua = v,
				.va = v
			}
		}
	);

	CHECK_EQ(res.size(), nx * ny * nz);

	for (size_t i = 0; i < nx * ny * nz; i++) {
		CHECK_EQ(res[i], nu);
	}
}
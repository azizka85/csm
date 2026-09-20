#include <doctest/doctest.h>

#include <generators/bathymetry/gebco_generator.h>

using namespace Generators::Bathymetry;

TEST_CASE("Test bathymetry GEBCO generator generate") {
	GEBCOGeneratorParams params{
		.bound = Data::GeoBound {
			.latMin = 36,
			.lonMin = 46,

			.latMax = 47,
			.lonMax = 56
		},
		.refDepth = -25,
		.minDepth = 1,
		.filePath = path(
			"data/bathymetry/gebco_2026_n47.0_s36.0_w46.0_e56.0.nc"
		)
	};

	GEBCOGenerator generator(params);

	auto depths = generator.generate(11, 11);

	CHECK_EQ(depths.size(), 121);
}
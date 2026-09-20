#include <doctest/doctest.h>

#include <generators/dz/triple_point_generator.h>

using namespace Calc::Grid;
using namespace Generators::DZ;

TEST_CASE("Test dz triple point generator generate") {
	double dzMin = 0.002;
	double dzMax = 0.1;		

	TriplePointGenerator generator(
		TriplePointGeneratorParams{
			.vertStep = make_unique<TriplePointVerticalStepSize>(
				TriplePointVerticalStepSizeParams{
					.zm = 0.51,
					.dzMin = dzMin,
					.dzMax = dzMax
				}
			),
			.dirAttr = ""
		}
	);

	auto res = generator.generate();

	CHECK_GT(res.size(), 0);

	CHECK_EQ(res[0], dzMin);

	for (size_t i = 0; i < res.size(); i++) {
		CHECK_LT(res[i], dzMin + dzMax);
	}
}
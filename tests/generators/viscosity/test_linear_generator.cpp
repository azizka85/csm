#include <doctest/doctest.h>

#include <generators/dz/uniform_generator.h>
#include <generators/bathymetry/uniform_generator.h>
#include <generators/viscosity/linear_generator.h>

using namespace Calc;
using namespace Generators;

TEST_CASE("Test viscosity linear generator params") {
	double h = 0;

	Turbulence::ViscosityState visc{
		.nut = 0,
		.nus = 0		
	};

	CHECK_THROWS_WITH_AS(
		Viscosity::LinearGenerator generator(h, visc),
		format("h should be > 0, but it is {}", h).c_str(),
		runtime_error
	);

	h = 100;

	CHECK_THROWS_WITH_AS(
		Viscosity::LinearGenerator generator(h, visc),
		format("nu should be > 0, but it is {}", visc.nus).c_str(),
		runtime_error
	);

	visc.nus = 0.4;

	CHECK_THROWS_WITH_AS(
		Viscosity::LinearGenerator generator(h, visc),
		format("nu should be > 0, but it is {}", visc.nut).c_str(),
		runtime_error
	);
}

TEST_CASE("Test viscosity linear generator generate") {
	size_t nx = 10;
	size_t ny = 11;

	double dz = 0.1;

	double h = 260;
	double ht = 100;

	Turbulence::ViscosityState visc{
		.nut = 0.005,
		.nus = 0.4
	};

	DZ::UniformGenerator dzGenerator(dz);
	Bathymetry::UniformGenerator hGenerator(h);
	Viscosity::LinearGenerator nuGenerator(ht, visc);

	auto dzData = dzGenerator.generate();

	size_t nz = dzData.size() + 1;

	auto bathymetry = hGenerator.generate(nx, ny, 1, 1);

	Data::VolumeState volume{
			.nx = nx,
			.ny = ny,
			.nz = nz,
			.dz = dzData,
			.bathymetry = bathymetry
	};

	vector<double> v;

	Data::WindSpeed windSpeed{
		.u10 = v,
		.v10 = v
	};

	Data::WindStress windStress{
		.qx = v,
		.qy = v
	};

	Data::Current current{
		.ua = v,
		.va = v
	};

	Data::WindCurrent windCurrent{
		.speed = windSpeed,
		.stress = windStress,
		.current = current
	};

	auto res = nuGenerator.generate(
		volume,
		windCurrent
	);

	CHECK_EQ(res.size(), nx * ny * nz);

	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
				auto p = j + i * ny;

				auto idt = p * nz;
				auto idb = nz - 1 + p * nz;

				CHECK_EQ(res[idt], visc.nus);
				CHECK_EQ(res[idb], visc.nut);
			}
	}

	ht = 300;

	nuGenerator = Viscosity::LinearGenerator(ht, visc);

	res = nuGenerator.generate(
		volume,
		windCurrent
	);

	CHECK_EQ(res.size(), nx * ny * nz);

	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			auto p = j + i * ny;

			auto idt = p * nz;
			auto idb = nz - 1 + p * nz;

			CHECK_EQ(res[idt], visc.nus);
			CHECK_GT(res[idb], visc.nut);
		}
	}
}
#include <doctest/doctest.h>

#include <generators/bathymetry/uniform_generator.h>
#include <generators/dz/uniform_generator.h>
#include <generators/wind/uniform_generator.h>
#include <generators/viscosity/wind_speed_generator.h>

using namespace Calc;
using namespace Generators;

TEST_CASE("Test viscosity wind speed generator generate") {
	Turbulence::Params turbulence{
		.k0 = 0.4,
		.k = 2e-5,
		.sigma = 1.2e-4,
		.num = 1.15e-6
	};

	Data::PhysicalParams phys{
		.rho = 1025,
		.f = 1.2e-4
	};

	auto viscosity = make_shared<Turbulence::DefaultViscosity>(
		Turbulence::DefaultViscosityParams{
			.turbulence = turbulence,
			.phys = phys
		}
	);

	Viscosity::WindSpeedGenerator nuGenerator(
		viscosity,
		viscosity,
		viscosity,
		""
	);

	size_t nx = 11;
	size_t ny = 11;

	double dz = 0.1;

	DZ::UniformGenerator dzGenerator(dz);

	auto dzData = dzGenerator.generate();

	size_t nz = dzData.size() + 1;

	double h = 260;

	Bathymetry::UniformGenerator hGenerator(h);

	auto hData = hGenerator.generate(nx, ny);

	double u10 = 22;
	double v10 = 0;

	double qx = 1.48225;
	double qy = 0;

	Generators::Wind::UniformGenerator qGenerator(
		Generators::Wind::UniformGeneratorParams{
			.speed = Calc::Wind::SpeedVector { .u10 = u10, .v10 = v10 },
			.stress = Calc::Wind::StressVector { .qx = qx, .qy = qy }
		}
	);

	auto windData = qGenerator.generate(nx, ny)[0];

	vector<double> ua(nx * ny, 0.5);
	vector<double> va(nx * ny, 0);

	Data::VolumeState volume{
			.nx = nx,
			.ny = ny,
			.nz = nz,
			.dx = 0,
			.dy = 0,
			.dz = dzData,
			.h = hData
	};

	Data::WindCurrent windCurrent{
		.speed = windData.speed,
		.stress = windData.stress,
		.current = Data::Current {
			.ua = ua,
			.va = va
		}
	};

	auto res = nuGenerator.generate(volume, windCurrent);

	CHECK_EQ(res.size(), nx * ny * nz);
}
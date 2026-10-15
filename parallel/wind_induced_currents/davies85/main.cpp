#include <iostream>

#include <calc/time.h>

#include <generators/area/proj_generator.h>
#include <generators/dz/triple_point_generator.h>
#include <generators/bathymetry/gebco_generator.h>
#include <generators/wind/ecmwf_generator.h>
#include <generators/viscosity/wind_speed_generator.h>

#include "double_precision/kernel.h"
#include "data.h"
#include "solver.h"

using namespace Parallel::WindInducedCurrents;

int main() {
	const double b = 1.1;

	const double f = 0.98e-4;
	const double g = 9.81;
	const double rho = 1025;
	const double kb = 0.002;

	const double k0 = 0.4;
	const double k = 2e-5;
	const double sigma = 1.2e-4;

	const double rhoAir = 1.225;
	const double Cd = 2.5e-3;

	const double latMin = 36;
	const double latMax = 47;

	const double lonMin = 46;
	const double lonMax = 56;

	double refDepth = -25;
	double minDepth = 5;

	const double dx = 13000;
	const double dy = 13000;

	const double dzMin = 0.002;
	const double dzMax = 0.1;

	const double zm = 0.51;

	const double ht = 100;

	const double num = 0.4;
	const double nut = 0.005;
	const double nu0 = 1.15e-6;

	const double w = 260;
	const double l = 260;

	const double hm = 260;

	const double u10m = 22;
	const double v10m = 22;

	const double qxm = 1.5;
	const double qym = 1.5;

	const double endTime = 90000;
	const double outputTimeStep = 600;

	const string gebcoFilePath = "data/bathymetry/gebco_2026_n47.0_s36.0_w46.0_e56.0.nc";
	const string ecmwfFilePath = "data/wind/ecmwf_2026_07_31_n47_e56_s36_w46.nc";

	const string outDir = "out";

	try {
		auto timeStep = make_shared<Calc::Time::DefaultAdjustTimeStep>(
			Calc::Time::DefaultAdjustTimeStepParams{
				.b = b,
				.tMax = outputTimeStep
			}
		);

		Calc::Data::GeoBound geoBound{
				.latMin = latMin,
				.lonMin = lonMin,

				.latMax = latMax,
				.lonMax = lonMax
		};

		auto areaGenerator = make_shared<Generators::Area::ProjGenerator>(geoBound);
		auto dzGenerator = make_shared<Generators::DZ::TriplePointGenerator>(
			Generators::DZ::TriplePointGeneratorParams{
				.vertStep = make_shared<Calc::Grid::TriplePointVerticalStepSize>(
					Calc::Grid::TriplePointVerticalStepSizeParams {
						.zm = zm,
						.dzMin = dzMin,
						.dzMax = dzMax
					}
				)
			}
		);

		auto hGenerator = make_shared<Generators::Bathymetry::GEBCOGenerator>(
			Generators::Bathymetry::GEBCOGeneratorParams{
				.bound = geoBound,
				.refDepth = refDepth,
				.minDepth = minDepth,
				.filePath = gebcoFilePath
			}
		);

		auto qGenerator = make_shared<Generators::Wind::ECMWFGenerator>(
			Generators::Wind::ECMWFGeneratorParams{
				.bound = geoBound,
				.windStress = make_shared<Calc::Wind::DefaultStress>(
					Calc::Wind::DefaultStressParams {
						.rhoAir = rhoAir,
						.cd = Cd
					}
				),
				.filePath = ecmwfFilePath
			}
		);

		PhysicalParams phys{ .rho = rho, .f = f };

		Calc::Turbulence::Params turbulence{
			.k0 = k0,
			.k = k,
			.sigma = sigma,
			.num = num
		};

		auto viscosity = make_shared<Calc::Turbulence::DefaultViscosity>(
			Calc::Turbulence::DefaultViscosityParams{
				.turbulence = turbulence,
				.phys = phys
			}
		);

		auto nuGenerator = make_shared<Generators::Viscosity::WindSpeedGenerator>(
			viscosity, viscosity, viscosity, ""
		);

		bool symmetrise = true;

		auto data = Davies85::Data(
			g, kb, k0, k, sigma, num,
			phys,
			dx, dy, symmetrise,
			areaGenerator,
			dzGenerator,
			hGenerator,
			qGenerator,
			nuGenerator
		);		

		data.init();

		vector<string> baseSources = {
			"opencl/double_precision/utils/max.cl",
			"opencl/double_precision/utils/interpolate.cl",
			"opencl/double_precision/slae/direct/tridiagonal.cl",
			"opencl/double_precision/wind_induced_currents/davies85/solver.cl",
		};

		vector<string> defaultViscositySources = {
			"opencl/double_precision/viscosity/default.cl"
		};

		vector<string> windSpeedViscositySources = {
			"opencl/double_precision/calc/turbulence.cl",
			"opencl/double_precision/calc/viscosity/wind_speed.cl",
			"opencl/double_precision/viscosity/wind_speed.cl"
		};

		auto sources = baseSources;

		sources.insert(sources.end(), windSpeedViscositySources.begin(), windSpeedViscositySources.end());

		bool calcNU = true;
		bool interpolateWind = true;

		auto kernel = make_shared<Davies85::DoublePrecision::Kernel>(
			data, 
			sources,
			interpolateWind
		);

		Davies85::Solver solver(
			data, 
			calcNU, interpolateWind,			
			endTime, outputTimeStep, outDir,
			timeStep,
			kernel
		);

		solver.solve();
	}
	catch (const exception& e) {
		cout << "Caught exception: " << e.what() << std::endl;
	}

	return 0;
}

#include <iostream>

#include <calc/time.h>

#include <generators/area/uniform_generator.h>
#include <generators/dz/triple_point_generator.h>
#include <generators/bathymetry/uniform_generator.h>
#include <generators/wind/uniform_generator.h>
#include <generators/viscosity/linear_generator.h>

#include "double_precision/kernel.h"
#include "data.h"
#include "solver.h"

using namespace Parallel::WindInducedCurrents;

int main() {
	const double b = 1.1;

	const double f = 1.2e-4;
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

	const double dx = 13;
	const double dy = 13;

	const double dzMin = 0.002;
	const double dzMax = 0.1;

	const double zm = 0.51;

	const double ht = 25;

	const double num = 0.4;
	const double nut = 0.05;
	const double nu0 = 1.15e-6;

	const double w = 260;
	const double l = 260;

	const double hm = 260;

	const double u10m = 0;
	const double v10m = -22;

	const double qxm = 0;
	const double qym = -1.5;

	const double tickTime = 600;
	const double endTime = 543600;
	const double outputTimeStep = 12600;

	const string gebcoFilePath = "data/bathymetry/gebco_2026_n45.5_s42.5_w48.0_e50.0.nc";
	const string ecmwfFilePath = "data/wind/ecmwf_2026_07_31_n45.5_e50_s42.5_w48.nc";

	const string outDir = "out";

	const bool symmetrise = false;

	const bool calcNU = false;
	const bool interpolateWind = false;

	bool writeOutput = true;

	try {
		auto timeStep = make_shared<Calc::Time::DefaultAdjustTimeStep>(
			Calc::Time::DefaultAdjustTimeStepParams{
				.b = b,
				.tMax = outputTimeStep
			}
		);		

		auto areaGenerator = make_shared<Generators::Area::UniformGenerator>(
			Generators::Area::Geometry{
				.l = l,
				.w = w
			}
		);
		auto dzGenerator = make_shared<Generators::DZ::TriplePointGenerator>(
			Generators::DZ::TriplePointGeneratorParams{
				.vertStep = make_shared<Calc::Grid::TriplePointVerticalStepSize>(
					Calc::Grid::TriplePointVerticalStepSizeParams{
						.zm = zm,
						.dzMin = dzMin,
						.dzMax = dzMax
					}
				)
			}
		);

		auto hGenerator = make_shared<Generators::Bathymetry::UniformGenerator>(hm);

		auto qGenerator = make_shared<Generators::Wind::UniformGenerator>(
			Generators::Wind::UniformGeneratorParams{
				.speed = Calc::Wind::SpeedVector{
					.u10 = u10m,
					.v10 = v10m
				},
				.stress = Calc::Wind::StressVector{
					.qx = qxm,
					.qy = qym
				}
			}
		);

		PhysicalParams phys{ .rho = rho, .f = f };
		
		auto nuGenerator = make_shared<Generators::Viscosity::LinearGenerator>(
			ht, 
			Calc::Turbulence::ViscosityState{
				.nut = nut,
				.nus = num
			}
		);		

		auto data = Davies85::Data(
			g, kb, k0, k, sigma, nu0,
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

		sources.insert(sources.end(), defaultViscositySources.begin(), defaultViscositySources.end());

		auto kernel = make_shared<Davies85::DoublePrecision::Kernel>(
			data, 
			sources,
			interpolateWind,
			1
		);		

		Davies85::Solver solver(
			data, 
			calcNU, interpolateWind, writeOutput,			
			tickTime, endTime, outputTimeStep, outDir,
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

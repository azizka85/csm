#ifndef PARALLEL_WIND_INDUCED_CURRENTS_DAVIES85_DATA_H
#define PARALLEL_WIND_INDUCED_CURRENTS_DAVIES85_DATA_H

#include <vector>
#include <memory>

#include <filesystem>

#include <calc/data.h>

#include <generators/area.h>
#include <generators/dz.h>
#include <generators/bathymetry.h>
#include <generators/wind.h>
#include <generators/viscosity.h>

using namespace std;
using namespace std::filesystem;

namespace Parallel::WindInducedCurrents::Davies85 {
	class Data {
	private:				
		shared_ptr<Generators::Area::IGenerator> areaGenerator;
		shared_ptr<Generators::DZ::IGenerator> dzGenerator;
		shared_ptr<Generators::Bathymetry::IGenerator> hGenerator;
		shared_ptr<Generators::Wind::IGenerator> qGenerator;
		shared_ptr<Generators::Viscosity::IGenerator> nuGenerator;

	public:
		size_t nx, ny, nz;
		double dx, dy;
		double rho, f, kb, g;
		double k0, k, sigma, num;

		vector<double> dz;
		vector<double> h;

		vector<double> ua;
		vector<double> va;

		vector<double> z;

		vector<double> uf;
		vector<double> vf;

		vector<double> ud;
		vector<double> vd;

		vector<double> nu;

		size_t currentWindIndex;
		vector<Calc::Data::Wind> windData;

		Data(
			double g, double kb, double k0, double k, double sigma, double num,
			Calc::Data::PhysicalParams phys,			
			double dx, double dy,
			shared_ptr<Generators::Area::IGenerator> areaGenerator,
			shared_ptr<Generators::DZ::IGenerator> dzGenerator,
			shared_ptr<Generators::Bathymetry::IGenerator> hGenerator,
			shared_ptr<Generators::Wind::IGenerator> qGenerator,
			shared_ptr<Generators::Viscosity::IGenerator> nuGenerator
		);

		void init();

		path addDirectory(path outDir);
	};
}

#endif
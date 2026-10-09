#ifndef PARALLEL_WIND_INDUCED_CURRENTS_DAVIES85_SOLVER_H
#define PARALLEL_WIND_INDUCED_CURRENTS_DAVIES85_SOLVER_H

#include <string>

#include <memory>

#include <calc/data.h>
#include <calc/time.h>

#include <generators/area.h>
#include <generators/dz.h>
#include <generators/bathymetry.h>
#include <generators/wind.h>
#include <generators/viscosity.h>

#include <writers/utils.h>
#include <writers/statistics.h>
#include <writers/surface.h>
#include <writers/volume.h>

using namespace std;

namespace Parallel::WindInducedCurrents::Davies85 {
	struct Directories {
		path root;
		path surface;
		path volume;
		path viscosity;
	};

	class Solver {
	private:
		PhysicalParams phys;

		double g;
		double kb;

		double dx;
		double dy;

		double endTime;
		double outputTimeStep;

		string outDir;

		shared_ptr<Calc::Time::IAdjustTimeStep> timeStep;
		shared_ptr<Generators::Area::IGenerator> areaGenerator;
		shared_ptr<Generators::DZ::IGenerator> dzGenerator;
		shared_ptr<Generators::Bathymetry::IGenerator> hGenerator;
		shared_ptr<Generators::Wind::IGenerator> qGenerator;
		shared_ptr<Generators::Viscosity::IGenerator> nuGenerator;

		Directories createDirectory();

		void setInitialCondition(
			int nx, int ny, int nz,
			vector<double>& uf,
			vector<double>& vf,
			vector<double>& ua,
			vector<double>& va,
			vector<double>& z
		);

		void writeData(
			const TimeState& time,
			const VolumeState& volume,
			Current& current,
			Calc::Data::Wind& currentWindData,
			vector<double>& z,
			const Directories& dirs,
			vector<double>& uf,
			vector<double>& ua,
			vector<double>& u,
			vector<double>& vf,
			vector<double>& va,
			vector<double>& v,
			vector<double>& w,
			vector<double>& nu,
			Writers::Surface::Writer& surfaceWriter,
			Writers::Volume::Writer& volumeElemWriter,
			Writers::Volume::Writer& volumeNodeWriter
		);

	public:
		Solver(
			PhysicalParams phys,
			double g, double kb,
			double dx, double dy,
			double endTime, double outputTimeStep, string outDir,
			shared_ptr<Calc::Time::IAdjustTimeStep> timeStep,
			shared_ptr<Generators::Area::IGenerator> areaGenerator,
			shared_ptr<Generators::DZ::IGenerator> dzGenerator,
			shared_ptr<Generators::Bathymetry::IGenerator> hGenerator,
			shared_ptr<Generators::Wind::IGenerator> qGenerator,
			shared_ptr<Generators::Viscosity::IGenerator> nuGenerator
		);

		void solve();
	};
}

#endif
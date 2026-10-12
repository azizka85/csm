#ifndef PARALLEL_WIND_INDUCED_CURRENTS_DAVIES85_SOLVER_H
#define PARALLEL_WIND_INDUCED_CURRENTS_DAVIES85_SOLVER_H

#include <calc/time.h>

#include <writers/utils.h>
#include <writers/statistics.h>
#include <writers/surface.h>
#include <writers/volume.h>

#include "data.h"
#include "kernel.h"

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
		Data& data;

		bool calcNU;

		double endTime;
		double outputTimeStep;

		string outDir;

		shared_ptr<Calc::Time::IAdjustTimeStep> timeStep;
		shared_ptr<IKernel> kernel;

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
			Calc::Data::WindCurrentData currentWindData,
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
			Data& data, bool calcNU,
			double endTime, double outputTimeStep, string outDir,
			shared_ptr<Calc::Time::IAdjustTimeStep> timeStep,
			shared_ptr<IKernel> kernel
		);

		void solve();
	};
}

#endif
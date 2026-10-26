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
		bool interpolateWind;

		bool writeOutput;

		double tickTime;
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

		void calculateSurfaceCurrent(
			int nx, int ny, int nz,
			const VolumeState& volume,
			const WindStress& windStress,
			const vector<double>& u,
			const vector<double>& v,
			const vector<double>& nu,
			vector<double>& us,
			vector<double>& vs
		);

		void calculateBottomStress(
			int nx, int ny, int nz,
			const VolumeState& volume,
			const vector<double>& uf, const vector<double>& vf,
			const vector<double>& ua, const vector<double>& va,
			const vector<double>& nu,
			vector<double>& tbx, vector<double>& tby
		);

		void writeData(
			const TimeState& time,
			const VolumeState& volume,
			const Calc::Data::WindCurrent& currentWindData,
			vector<double>& z,
			const Directories& dirs,
			vector<double>& uf,
			vector<double>& ua,
			vector<double>& us,
			vector<double>& u,
			vector<double>& vf,
			vector<double>& va,
			vector<double>& vs,
			vector<double>& v,
			vector<double>& w,
			vector<double>& nu,
			vector<double>& tbx,
			vector<double>& tby,
			Writers::Surface::Writer& surfaceWriter,
			Writers::Volume::Writer& volumeElemWriter,
			Writers::Volume::Writer& volumeNodeWriter
		);

	public:
		Solver(
			Data& data, 
			bool calcNU, bool interpolateWind, bool writeOutput,
			double tickTime, double endTime, double outputTimeStep, string outDir,
			shared_ptr<Calc::Time::IAdjustTimeStep> timeStep,
			shared_ptr<IKernel> kernel
		);

		void solve();
	};
}

#endif
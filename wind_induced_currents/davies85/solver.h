#ifndef WIND_INDUCED_CURRENTS_DAVIES85_SOLVER_H
#define WIND_INDUCED_CURRENTS_DAVIES85_SOLVER_H

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

namespace WindInducedCurrents::Davies85 {
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

		void calcUA(
			int nx, int ny, int nz, 
			PhysicalParams phys, 
			double kb, double g, 
			double dx, double dt, 
			vector<double>& h,
			vector<double>& qx,
			vector<double>& ua,
			vector<double>& va,
			vector<double>& z,
			vector<double>& uf,
			vector<double>& u1a
		);

		void calcVA(
			int nx, int ny, int nz,
			PhysicalParams phys,
			double kb, double g,
			double dy, double dt,
			vector<double>& h,
			vector<double>& qy,
			vector<double>& ua,
			vector<double>& va,
			vector<double>& z,
			vector<double>& vf,
			vector<double>& v1a
		);

		void calcRHS(
			int nx, int ny, int nz,
			PhysicalParams phys, double kb, 
			double dx, double dy, double dt,
			vector<double>& dz, 
			vector<double>& h,
			vector<double>& nu,
			vector<double>& qx,
			vector<double>& qy,
			vector<double>& ua,
			vector<double>& u1a,
			vector<double>& va,
			vector<double>& v1a,
			vector<double>& uf,
			vector<double>& vf,
			vector<double>& ud,
			vector<double>& vd
		);

		void createTridiagonalMatrix(
			double kb, double dt,
			int nx, int ny, int nz,
			vector<double>& dz,
			vector<double>& h,
			vector<double>& nu,
			vector<double>& al,
			vector<double>& ac,
			vector<double>& ar
		);

		void calcUVF(
			int nx, int ny, int nz,
			vector<double>& h,
			vector<double>& al,
			vector<double>& ac,
			vector<double>& ar,
			vector<double>& uf,
			vector<double>& vf,
			vector<double>& ud,
			vector<double>& vd
		);

		void calcZ(
			int nx, int ny, int nz,
			double dt, double dx, double dy,
			vector<double>& h,
			vector<double>& ua,
			vector<double>& va,
			vector<double>& z
		);

		double maxAbsResidual(
			int nx, int ny, int nz,
			vector<double>& h,
			vector<double>& al,
			vector<double>& ac,
			vector<double>& ar,
			vector<double>& uf,
			vector<double>& vf,
			vector<double>& ud,
			vector<double>& vd
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

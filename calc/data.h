#ifndef CALC_DATA_H
#define CALC_DATA_H

#include <string>
#include <vector>

using namespace std;

namespace Calc::Data {
	struct PhysicalParams {
		double rho;
		double f;		
	};

	struct DepthAveragedVelocity {
		double ua;
		double va;
	};

	struct Column {
		string name;
		size_t precision;
	};

	struct TimeState {
		size_t m;
		double t;
	};

	struct SurfaceState {
		size_t nx;
		size_t ny;

		double dx;
		double dy;
	};

	struct SurfaceVector {
		vector<double> &vx;
		vector<double> &vy;
	};

	struct VolumeState {
		size_t nx;
		size_t ny;
		size_t nz;

		double dx;
		double dy;

		vector<double> &dz;
		vector<double> &h;
	};

	struct VolumeVector {
		vector<double> &vx;
		vector<double> &vy;
		vector<double> &vz;
	};

	struct WindSpeed {
		vector<double> u10;
		vector<double> v10;
	};

	struct WindStress {
		vector<double> qx;
		vector<double> qy;
	};

	struct Wind {
		double time;
		WindSpeed speed;
		WindStress stress;
	};

	struct Current {
		vector<double> &ua;
		vector<double> &va;
	};

	struct WindCurrent {
		WindSpeed speed;
		WindStress stress;
		Current current;
	};

	struct CartesianBound {
		double minX;
		double minY;

		double maxX;
		double maxY;
	};

	struct GeoBound {
		double latMin;
		double lonMin;

		double latMax;
		double lonMax;
	};

	struct CartesianPoint {
		double x;
		double y;
	};

	struct GeoPoint {
		double lat;
		double lon;
	};
}

#endif 
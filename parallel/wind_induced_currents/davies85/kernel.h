#ifndef PARALLEL_WIND_INDUCED_CURRENTS_DAVIES85_KERNEL_H
#define PARALLEL_WIND_INDUCED_CURRENTS_DAVIES85_KERNEL_H

#include <filesystem>

using namespace std::filesystem;

namespace Parallel::WindInducedCurrents::Davies85 {
	class IKernel {
	public:
		virtual void calcUA(double dt) = 0;
		virtual void calcVA(double dt) = 0;

		virtual void updateUVA() = 0;
		virtual void swapUA() = 0;
		virtual void swapVA() = 0;

		virtual void calcRHSU(double dt) = 0;
		virtual void calcUF() = 0;

		virtual void calcRHSV(double dt) = 0;
		virtual void calcVF() = 0;

		virtual void createTridiagonalMatrix(double dt) = 0;
		virtual void calcResidual() = 0;		

		virtual void updateUVF() = 0;
		virtual void updateUVD() = 0;

		virtual void calcZ(double dt) = 0;
		virtual void updateZ() = 0;

		virtual void updateQUV10Device(bool interpolateWind) = 0;
		virtual void interpolateQUV10(double alpha) = 0;
		
		virtual void updateQUV10(
			vector<double>& u10, vector<double>& v10,
			vector<double>& qx, vector<double>& qy
		) = 0;

		virtual void calcNU() = 0;
		virtual double calcNUMax() = 0;
		virtual void updateNU() = 0;

		virtual path addDirectory(path outDir) = 0;
	};
}

#endif
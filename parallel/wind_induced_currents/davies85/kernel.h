#ifndef PARALLEL_WIND_INDUCED_CURRENTS_DAVIES85_KERNEL_H
#define PARALLEL_WIND_INDUCED_CURRENTS_DAVIES85_KERNEL_H

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include <vector>

using namespace std;

namespace Parallel {
	namespace WindInducedCurrents {
		namespace Davies85 {
			class Kernel {
			private:
				dim3 threadsPerBlock;
				dim3 blocksPerGrid2D;
				dim3 blocksPerGrid3D;

				size_t nx;
				size_t ny;
				size_t nz;

				double dx;
				double dy;

				double rho;
				double f;

				double kb; 
				double g;

				vector<double>& dz;
				vector<double>& h;

				vector<double>& ua;
				vector<double>& va;

				vector<double>& z;

				vector<double>& uf;
				vector<double>& vf;

				vector<double>& u10;
				vector<double>& v10;

				vector<double>& qx;
				vector<double>& qy;

				vector<double>& nu;

				double* d_dz;
				double* d_h;

				double* d_ua;
				double* d_u1a;

				double* d_va;
				double* d_v1a;

				double* d_z;
				double* d_zt;

				double* d_uf;
				double* d_uft;
				double* d_ud;
				double* d_udt;

				double* d_vf;
				double* d_vft;
				double* d_vd;
				double* d_vdt;

				double* d_al;
				double* d_ac;
				double* d_ar;

				double* d_u10;
				double* d_v10;

				double* d_qx;
				double* d_qy;

				double* d_nu;

				double* d_alpha;

				double* d_res;
				double* d_res2;

				int* d_res_idx;

			public:
				Kernel(
					size_t nx, size_t ny,
					double dx, double dy, 
					double rho, double f, double kb, double g, 
					vector<double>& dz, vector<double>& h,
					vector<double>& ua, vector<double>& va,
					vector<double>& z,
					vector<double>& uf, vector<double>& vf,
					vector<double>& u10, vector<double>& v10,
					vector<double>& qx, vector<double>& qy,
					vector<double>& nu
				);

				void calcUA(double dt);
				void calcVA(double dt);
				
				void updateUVA();
				void swapUA();
				void swapVA();

				void calcRHS(double dt);
				void createTridiagonalMatrix(double dt);
				
				void calcUVF();
				void updateUVF();
				void updateUVFTDevice();
				double maxAbsUFDifference();
				double maxAbsVFDifference();

				void calcZ(double dt);
				void updateZ();
				void updateZTDevice();
				double maxAbsZDifference();

				double maxAbsResidual();

				void updateQUV10Device();

				void calcNU();
				void updateNU();
			};
		}
	}
}

#endif
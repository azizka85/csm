#ifndef PARALLEL_WIND_INDUCED_CURRENTS_DAVIES85_SINGLE_PRECISION_KERNEL_H
#define PARALLEL_WIND_INDUCED_CURRENTS_DAVIES85_SINGLE_PRECISION_KERNEL_H

#include <CL/opencl.hpp>

#include <vector>

#include "../data.h"
#include "../kernel.h"

using namespace std;
using namespace std::filesystem;

namespace Parallel::WindInducedCurrents::Davies85::SinglePrecision {
	class Kernel: public IKernel {
	private:
		cl::NDRange surfaceRange;
		cl::NDRange volumeRange;

		cl::NDRange valueRange;
		cl::NDRange rowRange;

		Davies85::Data& data;

		cl::Buffer bufferDZ;
		cl::Buffer bufferH;

		cl::Buffer bufferUA;
		cl::Buffer bufferU1A;

		cl::Buffer bufferVA;
		cl::Buffer bufferV1A;

		cl::Buffer bufferZ;

		cl::Buffer bufferUF;
		cl::Buffer bufferUD;
		cl::Buffer bufferUDT;

		cl::Buffer bufferVF;
		cl::Buffer bufferVD;
		cl::Buffer bufferVDT;

		cl::Buffer bufferAL;
		cl::Buffer bufferAC;
		cl::Buffer bufferAR;

		cl::Buffer bufferU10;
		cl::Buffer bufferV10;

		cl::Buffer bufferQX;
		cl::Buffer bufferQY;

		cl::Buffer bufferNU;

		cl::Buffer bufferV;
		cl::Buffer bufferR;
		cl::Buffer bufferP;
						
		cl::CommandQueue queue;

		cl::Kernel calcUAKernel;
		cl::Kernel calcVAKernel;

		cl::Kernel calcRHSKernel;
		cl::Kernel createTridiagonalMatrixKernel;
		cl::Kernel calcUVFKernel;
		cl::Kernel calcResidualKernel;

		cl::Kernel calcZKernel;

		cl::Kernel calcNUKernel;

		cl::Kernel calcMaxPlaneKernel;
		cl::Kernel calcMaxRowKernel;
		cl::Kernel calcMaxKernel;

		vector<float> scratch;

		string loadKernelSource(path filePath);
		void loadKernelSources(cl::Program::Sources& sources, vector<string>& sourcePaths);

		void readFloatBuffer(const cl::Buffer& buf, vector<double>& dst, size_t count);

		void printDeviceDescription(cl::Platform& platform, cl::Device& device);

	public:
		Kernel(
			Davies85::Data& data,
			vector<string>& sourcePaths,
			size_t platformIdx = 0, size_t deviceIdx = 0
		);

		void calcUA(double dt) override;
		void calcVA(double dt) override;

		void updateUVA() override;
		void swapUA() override;
		void swapVA() override;

		void calcRHS(double dt) override;
		void createTridiagonalMatrix(double dt) override;
		void calcResidual() override;

		void calcUVF() override;
		void updateUVF() override;
		void updateUVD() override;

		void calcZ(double dt) override;
		void updateZ() override;

		void updateQUV10Device() override;

		void calcNU() override;
		double calcNUMax() override;
		void updateNU() override;

		path addDirectory(path outDir) override;		
};
}

#endif
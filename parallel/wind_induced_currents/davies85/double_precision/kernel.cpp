#include <fstream>

#include <stdexcept>

#include "kernel.h"

using namespace Parallel::WindInducedCurrents::Davies85::DoublePrecision;

Kernel::Kernel(
	Davies85::Data& data,
	string viscosityKernelSource,
	string viscosityFuncSource
): data(data)
{
	surfaceRange = cl::NDRange(data.nx, data.ny);
	volumeRange = cl::NDRange(data.nx, data.ny, data.nz);

	vector<cl::Platform> platforms;

	cl::Platform::get(&platforms);

	cl::Platform platform = platforms.front();

	vector<cl::Device> devices;

	platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);

	cl::Device device = devices.front();
	cl::Context context(device);
	queue = cl::CommandQueue(context, device);

	cl::Program::Sources sources;

	loadKernelSources(sources, viscosityKernelSource, viscosityFuncSource);

	cl::Program program(context, sources);

	auto err = program.build(device);

	if (err == -11) {
		auto build_log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>();

		bool ok = true;
	}

	bufferDZ = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nz);
	err = queue.enqueueWriteBuffer(bufferDZ, CL_TRUE, 0, sizeof(double) * data.nz, data.dz.data());

	bufferH = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferH, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.h.data());

	bufferUA = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferUA, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.ua.data());

	bufferU1A = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);

	bufferVA = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferVA, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.va.data());

	bufferV1A = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);

	bufferZ = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferZ, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.z.data());

	bufferUF = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny * data.nz);
	err = queue.enqueueWriteBuffer(bufferUF, CL_TRUE, 0, sizeof(double) * data.nx * data.ny * data.nz, data.uf.data());

	bufferUD = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny * data.nz);
	bufferUDT = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny * data.nz);

	bufferVF = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny * data.nz);
	err = queue.enqueueWriteBuffer(bufferVF, CL_TRUE, 0, sizeof(double) * data.nx * data.ny * data.nz, data.vf.data());

	bufferVD = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny * data.nz);
	bufferVDT = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny * data.nz);

	bufferAL = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny * data.nz);
	bufferAC = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny * data.nz);
	bufferAR = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny * data.nz);

	bufferU10 = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferU10, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex].speed.u10.data());

	bufferV10 = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferV10, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex].speed.v10.data());

	bufferQX = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferQX, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex].stress.qx.data());

	bufferQY = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferQY, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex].stress.qy.data());

	bufferNU = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny * (data.nz + 1));
	err = queue.enqueueWriteBuffer(bufferNU, CL_TRUE, 0, sizeof(double) * data.nx * data.ny * (data.nz + 1), data.nu.data());

	calcUAKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_ua");

	calcUAKernel.setArg(0, data.f);
	calcUAKernel.setArg(1, data.kb);
	calcUAKernel.setArg(2, data.g);
	calcUAKernel.setArg(3, data.rho);
	calcUAKernel.setArg(4, data.dx);
	calcUAKernel.setArg(6, data.nx);
	calcUAKernel.setArg(7, data.ny);
	calcUAKernel.setArg(8, data.nz);
	calcUAKernel.setArg(9, bufferH);
	calcUAKernel.setArg(10, bufferQX);
	calcUAKernel.setArg(11, bufferUA);
	calcUAKernel.setArg(12, bufferVA);
	calcUAKernel.setArg(13, bufferZ);
	calcUAKernel.setArg(14, bufferUF);
	calcUAKernel.setArg(15, bufferU1A);

	calcVAKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_va");

	calcVAKernel.setArg(0, data.f);
	calcVAKernel.setArg(1, data.kb);
	calcVAKernel.setArg(2, data.g);
	calcVAKernel.setArg(3, data.rho);
	calcVAKernel.setArg(4, data.dy);
	calcVAKernel.setArg(6, data.nx);
	calcVAKernel.setArg(7, data.ny);
	calcVAKernel.setArg(8, data.nz);
	calcVAKernel.setArg(9, bufferH);
	calcVAKernel.setArg(10, bufferQY);
	calcVAKernel.setArg(11, bufferUA);
	calcVAKernel.setArg(12, bufferVA);
	calcVAKernel.setArg(13, bufferZ);
	calcVAKernel.setArg(14, bufferVF);
	calcVAKernel.setArg(15, bufferV1A);

	calcRHSKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_rhs");

	calcRHSKernel.setArg(0, data.f);
	calcRHSKernel.setArg(1, data.kb);
	calcRHSKernel.setArg(2, data.rho);
	calcRHSKernel.setArg(3, data.dx);
	calcRHSKernel.setArg(4, data.dy);
	calcRHSKernel.setArg(6, data.ny);
	calcRHSKernel.setArg(7, data.nz);
	calcRHSKernel.setArg(8, bufferNU);
	calcRHSKernel.setArg(9, bufferDZ);
	calcRHSKernel.setArg(10, bufferH);
	calcRHSKernel.setArg(11, bufferQX);
	calcRHSKernel.setArg(12, bufferQY);
	calcRHSKernel.setArg(13, bufferUA);
	calcRHSKernel.setArg(14, bufferVA);
	calcRHSKernel.setArg(15, bufferUF);
	calcRHSKernel.setArg(16, bufferVF);
	calcRHSKernel.setArg(17, bufferUD);
	calcRHSKernel.setArg(18, bufferVD);

	createTridiagonalMatrixKernel = cl::Kernel(program, "wind_induced_currents_davies85_create_tridiagonal_matrix");

	createTridiagonalMatrixKernel.setArg(0, data.kb);
	createTridiagonalMatrixKernel.setArg(2, data.ny);
	createTridiagonalMatrixKernel.setArg(3, data.nz);
	createTridiagonalMatrixKernel.setArg(4, bufferNU);
	createTridiagonalMatrixKernel.setArg(5, bufferDZ);
	createTridiagonalMatrixKernel.setArg(6, bufferH);
	createTridiagonalMatrixKernel.setArg(7, bufferAL);
	createTridiagonalMatrixKernel.setArg(8, bufferAC);
	createTridiagonalMatrixKernel.setArg(9, bufferAR);

	calcUVFKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_uvf");

	calcUVFKernel.setArg(0, data.ny);
	calcUVFKernel.setArg(1, data.nz);
	calcUVFKernel.setArg(2, bufferH);
	calcUVFKernel.setArg(3, bufferAL);
	calcUVFKernel.setArg(4, bufferAC);
	calcUVFKernel.setArg(5, bufferAR);
	calcUVFKernel.setArg(6, bufferUF);
	calcUVFKernel.setArg(7, bufferVF);
	calcUVFKernel.setArg(8, bufferUD);
	calcUVFKernel.setArg(9, bufferVD);

	calcResidualKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_residual");

	calcResidualKernel.setArg(0, data.ny);
	calcResidualKernel.setArg(1, data.nz);
	calcResidualKernel.setArg(2, bufferH);
	calcResidualKernel.setArg(3, bufferAL);
	calcResidualKernel.setArg(4, bufferAC);
	calcResidualKernel.setArg(5, bufferAR);
	calcResidualKernel.setArg(6, bufferUF);
	calcResidualKernel.setArg(7, bufferVF);
	calcResidualKernel.setArg(8, bufferUDT);
	calcResidualKernel.setArg(9, bufferVDT);

	calcZKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_z");

	calcZKernel.setArg(0, data.dx);
	calcZKernel.setArg(1, data.dy);
	calcZKernel.setArg(3, data.nx);
	calcZKernel.setArg(4, data.ny);
	calcZKernel.setArg(5, data.nz);
	calcZKernel.setArg(6, bufferH);
	calcZKernel.setArg(7, bufferUA);
	calcZKernel.setArg(8, bufferVA);
	calcZKernel.setArg(9, bufferZ);
}

string Kernel::loadKernelSource(path filePath) {
	ifstream file(filePath);

	if (file.bad()) throw runtime_error("Cannot open file: " + filePath.string());

	return string(
		istreambuf_iterator<char>(file),
		istreambuf_iterator<char>()
	);
}

void Kernel::loadKernelSources(
	cl::Program::Sources& sources,
	string viscosityKernelSource,
	string viscosityFuncSource
) {
	sources.push_back(
		loadKernelSource(
			"opencl/double_precision/wind_induced_currents/davies85/solver.cl"
		)
	);

	sources.push_back(
		loadKernelSource(
			viscosityFuncSource
		)
	);

	sources.push_back(
		loadKernelSource(
			viscosityKernelSource
		)
	);
}

void Kernel::calcUA(double dt) {
	calcUAKernel.setArg(5, dt);
	calcUAKernel.setArg(11, bufferUA);
	calcUAKernel.setArg(12, bufferVA);
	calcUAKernel.setArg(15, bufferU1A);

	auto err = queue.enqueueNDRangeKernel(calcUAKernel, cl::NullRange, surfaceRange, cl::NullRange);	
}

void Kernel::calcVA(double dt) {
	calcVAKernel.setArg(5, dt);
	calcVAKernel.setArg(11, bufferUA);
	calcVAKernel.setArg(12, bufferVA);
	calcVAKernel.setArg(15, bufferV1A);

	auto err = queue.enqueueNDRangeKernel(calcVAKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

void Kernel::updateUVA() {
	auto err = queue.enqueueReadBuffer(bufferUA, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.ua.data());
	err = queue.enqueueReadBuffer(bufferVA, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.va.data());
}

void Kernel::swapUA() {
	swap(bufferU1A, bufferUA);
}

void Kernel::swapVA() {
	swap(bufferV1A, bufferVA);
}

void Kernel::calcRHS(double dt) {
	calcRHSKernel.setArg(5, dt);
	calcRHSKernel.setArg(13, bufferUA);
	calcRHSKernel.setArg(14, bufferVA);

	auto err = queue.enqueueNDRangeKernel(calcRHSKernel, cl::NullRange, volumeRange, cl::NullRange);

	err = queue.enqueueCopyBuffer(bufferUD, bufferUDT, 0, 0, sizeof(double) * data.nx * data.ny * data.nz);
	err = queue.enqueueCopyBuffer(bufferVD, bufferVDT, 0, 0, sizeof(double) * data.nx * data.ny * data.nz);
}

void Kernel::createTridiagonalMatrix(double dt) {
	createTridiagonalMatrixKernel.setArg(1, dt);

	auto err = queue.enqueueNDRangeKernel(createTridiagonalMatrixKernel, cl::NullRange, volumeRange, cl::NullRange);
}

void Kernel::calcResidual() {
	auto err = queue.enqueueNDRangeKernel(calcResidualKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

void Kernel::calcUVF() {
	auto err = queue.enqueueNDRangeKernel(calcUVFKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

void Kernel::updateUVF() {
	auto err = queue.enqueueReadBuffer(bufferUF, CL_TRUE, 0, sizeof(double) * data.nx * data.ny * data.nz, data.uf.data());
	err = queue.enqueueReadBuffer(bufferVF, CL_TRUE, 0, sizeof(double) * data.nx * data.ny * data.nz, data.vf.data());
}

void Kernel::updateUVD() {
	auto err = queue.enqueueReadBuffer(bufferUDT, CL_TRUE, 0, sizeof(double) * data.nx * data.ny * data.nz, data.ud.data());
	err = queue.enqueueReadBuffer(bufferVDT, CL_TRUE, 0, sizeof(double) * data.nx * data.ny * data.nz, data.vd.data());
}

void Kernel::calcZ(double dt) {
	calcZKernel.setArg(2, dt);
	calcZKernel.setArg(7, bufferUA);
	calcZKernel.setArg(8, bufferVA);

	auto err = queue.enqueueNDRangeKernel(calcZKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

void Kernel::updateZ() {
	auto err = queue.enqueueReadBuffer(bufferZ, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.z.data());
}

void Kernel::updateQUV10Device() {
	auto err = queue.enqueueWriteBuffer(bufferU10, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex].speed.u10.data());
	err = queue.enqueueWriteBuffer(bufferV10, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex].speed.v10.data());

	err = queue.enqueueWriteBuffer(bufferQX, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex].stress.qx.data());
	err = queue.enqueueWriteBuffer(bufferQY, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex].stress.qy.data());
}

void Kernel::calcNU() {}

void Kernel::updateNU() {}

path Kernel::addDirectory(path outDir) {
	return outDir / path("dp");
}

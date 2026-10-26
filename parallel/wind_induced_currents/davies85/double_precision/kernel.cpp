#include <iostream>
#include <fstream>

#include <stdexcept>

#include "kernel.h"

using namespace Parallel::WindInducedCurrents::Davies85::DoublePrecision;

Kernel::Kernel(
	Davies85::Data& data,
	vector<string>& sourcePaths,
	bool interpolateWind,
	size_t platformIdx, size_t deviceIdx
): data(data)
{
	surfaceRange = cl::NDRange(data.nx, data.ny);
	volumeRange = cl::NDRange(data.nx, data.ny, data.nz);

	valueRange = cl::NDRange(1);
	rowRange = cl::NDRange(data.nx);

	vector<cl::Platform> platforms;

	cl::Platform::get(&platforms);

	cl::Platform platform = platforms[platformIdx];

	vector<cl::Device> devices;

	platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);

	cl::Device device = devices[deviceIdx];

	printDeviceDescription(platform, device);

	cl::Context context(device);

	queue = cl::CommandQueue(context, device);

	cl::Program::Sources sources;

	loadKernelSources(sources, sourcePaths);

	cl::Program program(context, sources);

	auto err = program.build(device);

	if (err == -11) {
		auto build_log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>();

		bool ok = true;
	}

	bufferDZ = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nz);
	err = queue.enqueueWriteBuffer(bufferDZ, CL_TRUE, 0, sizeof(double) * data.nz, data.dz.data());

	bufferH = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferH, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.bathymetry.h.data());

	bufferUA = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferUA, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.ua.data());

	bufferU1A = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueCopyBuffer(bufferUA, bufferU1A, 0, 0, sizeof(double) * data.nx * data.ny);	

	bufferVA = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferVA, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.va.data());

	bufferV1A = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
	err = queue.enqueueCopyBuffer(bufferVA, bufferV1A, 0, 0, sizeof(double) * data.nx * data.ny);

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

	if (interpolateWind) {
		bufferU10C = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
		err = queue.enqueueCopyBuffer(bufferU10, bufferU10C, 0, 0, sizeof(double) * data.nx * data.ny);

		bufferV10C = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
		err = queue.enqueueCopyBuffer(bufferV10, bufferV10C, 0, 0, sizeof(double) * data.nx * data.ny);

		bufferQXC = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
		err = queue.enqueueCopyBuffer(bufferQX, bufferQXC, 0, 0, sizeof(double) * data.nx * data.ny);

		bufferQYC = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
		err = queue.enqueueCopyBuffer(bufferQY, bufferQYC, 0, 0, sizeof(double) * data.nx * data.ny);

		bufferU10N = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
		bufferV10N = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
		bufferQXN = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);
		bufferQYN = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);

		if (data.currentWindIndex + 1 < data.windData.size()) {
			err = queue.enqueueWriteBuffer(bufferU10N, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex + 1].speed.u10.data());
			err = queue.enqueueWriteBuffer(bufferV10N, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex + 1].speed.v10.data());
			err = queue.enqueueWriteBuffer(bufferQXN, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex + 1].stress.qx.data());
			err = queue.enqueueWriteBuffer(bufferQYN, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex + 1].stress.qy.data());
		}
		else {
			err = queue.enqueueCopyBuffer(bufferU10, bufferU10N, 0, 0, sizeof(double) * data.nx * data.ny);
			err = queue.enqueueCopyBuffer(bufferV10, bufferV10N, 0, 0, sizeof(double) * data.nx * data.ny);
			err = queue.enqueueCopyBuffer(bufferQX, bufferQXN, 0, 0, sizeof(double) * data.nx * data.ny);
			err = queue.enqueueCopyBuffer(bufferQY, bufferQYN, 0, 0, sizeof(double) * data.nx * data.ny);
		}
	}

	bufferNU = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny * (data.nz + 1));
	err = queue.enqueueWriteBuffer(bufferNU, CL_TRUE, 0, sizeof(double) * data.nx * data.ny * (data.nz + 1), data.nu.data());

	bufferV = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double));
	bufferR = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx);
	bufferP = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double) * data.nx * data.ny);

	calcUAKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_ua");

	calcUAKernel.setArg(0, data.f);
	calcUAKernel.setArg(1, data.kb);
	calcUAKernel.setArg(2, data.g);
	calcUAKernel.setArg(3, data.rho);
	calcUAKernel.setArg(4, data.dx);
	calcUAKernel.setArg(6, (int)data.nx);
	calcUAKernel.setArg(7, (int)data.ny);
	calcUAKernel.setArg(8, (int)data.nz);
	calcUAKernel.setArg(9, bufferH);
	calcUAKernel.setArg(10, bufferDZ);
	calcUAKernel.setArg(11, bufferQX);
	calcUAKernel.setArg(12, bufferUA);
	calcUAKernel.setArg(13, bufferVA);
	calcUAKernel.setArg(14, bufferZ);
	calcUAKernel.setArg(15, bufferUF);
	calcUAKernel.setArg(16, bufferNU);
	calcUAKernel.setArg(17, bufferU1A);

	calcVAKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_va");

	calcVAKernel.setArg(0, data.f);
	calcVAKernel.setArg(1, data.kb);
	calcVAKernel.setArg(2, data.g);
	calcVAKernel.setArg(3, data.rho);
	calcVAKernel.setArg(4, data.dy);
	calcVAKernel.setArg(6, (int)data.nx);
	calcVAKernel.setArg(7, (int)data.ny);
	calcVAKernel.setArg(8, (int)data.nz);
	calcVAKernel.setArg(9, bufferH);
	calcVAKernel.setArg(10, bufferDZ);
	calcVAKernel.setArg(11, bufferQY);
	calcVAKernel.setArg(12, bufferUA);
	calcVAKernel.setArg(13, bufferVA);
	calcVAKernel.setArg(14, bufferZ);
	calcVAKernel.setArg(15, bufferVF);
	calcVAKernel.setArg(16, bufferNU);
	calcVAKernel.setArg(17, bufferV1A);	

	createTridiagonalMatrixKernel = cl::Kernel(program, "wind_induced_currents_davies85_create_tridiagonal_matrix");

	createTridiagonalMatrixKernel.setArg(0, data.kb);
	createTridiagonalMatrixKernel.setArg(2, (int)data.ny);
	createTridiagonalMatrixKernel.setArg(3, (int)data.nz);
	createTridiagonalMatrixKernel.setArg(4, bufferNU);
	createTridiagonalMatrixKernel.setArg(5, bufferDZ);
	createTridiagonalMatrixKernel.setArg(6, bufferH);
	createTridiagonalMatrixKernel.setArg(7, bufferAL);
	createTridiagonalMatrixKernel.setArg(8, bufferAC);
	createTridiagonalMatrixKernel.setArg(9, bufferAR);

	calcRHSUKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_rhs_u");
	calcRHSUKernel.setArg(0, data.f);
	calcRHSUKernel.setArg(1, data.kb);
	calcRHSUKernel.setArg(2, data.rho);
	calcRHSUKernel.setArg(4, (int)data.ny);
	calcRHSUKernel.setArg(5, (int)data.nz);
	calcRHSUKernel.setArg(6, bufferNU);
	calcRHSUKernel.setArg(7, bufferDZ);
	calcRHSUKernel.setArg(8, bufferH);
	calcRHSUKernel.setArg(9, bufferQX);
	calcRHSUKernel.setArg(11, bufferUF);
	calcRHSUKernel.setArg(12, bufferVF);
	calcRHSUKernel.setArg(13, bufferUD);

	calcRHSVKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_rhs_v");
	calcRHSVKernel.setArg(0, data.f);
	calcRHSVKernel.setArg(1, data.kb);
	calcRHSVKernel.setArg(2, data.rho);
	calcRHSVKernel.setArg(4, (int)data.ny);
	calcRHSVKernel.setArg(5, (int)data.nz);
	calcRHSVKernel.setArg(6, bufferNU);
	calcRHSVKernel.setArg(7, bufferDZ);
	calcRHSVKernel.setArg(8, bufferH);
	calcRHSVKernel.setArg(9, bufferQY);
	calcRHSVKernel.setArg(11, bufferUF);
	calcRHSVKernel.setArg(12, bufferVF);	
	calcRHSVKernel.setArg(13, bufferVD);

	calcUFKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_column");
	calcVFKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_column");

	for (auto* kr : { &calcUFKernel, &calcVFKernel }) {
		kr->setArg(0, (int)data.ny);
		kr->setArg(1, (int)data.nz);
		kr->setArg(2, bufferH);
		kr->setArg(3, bufferDZ);
		kr->setArg(4, bufferAL);
		kr->setArg(5, bufferAC);
		kr->setArg(6, bufferAR);
	}

	calcUFKernel.setArg(7, bufferUD);  
	calcUFKernel.setArg(8, bufferUF);

	calcVFKernel.setArg(7, bufferVD);  
	calcVFKernel.setArg(8, bufferVF);

	calcResidualKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_residual");

	calcResidualKernel.setArg(0, (int)data.ny);
	calcResidualKernel.setArg(1, (int)data.nz);
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
	calcZKernel.setArg(3, (int)data.nx);
	calcZKernel.setArg(4, (int)data.ny);
	calcZKernel.setArg(5, (int)data.nz);
	calcZKernel.setArg(6, bufferH);
	calcZKernel.setArg(7, bufferUA);
	calcZKernel.setArg(8, bufferVA);
	calcZKernel.setArg(9, bufferZ);

	calcNUKernel = cl::Kernel(program, "viscosity");

	calcNUKernel.setArg(0, data.k0);
	calcNUKernel.setArg(1, data.rho);
	calcNUKernel.setArg(2, data.f);
	calcNUKernel.setArg(3, data.k);
	calcNUKernel.setArg(4, data.sigma);
	calcNUKernel.setArg(5, data.num);
	calcNUKernel.setArg(6, (int)data.ny);
	calcNUKernel.setArg(7, (int)data.nz);
	calcNUKernel.setArg(8, bufferDZ);
	calcNUKernel.setArg(9, bufferH);
	calcNUKernel.setArg(10, bufferU10);
	calcNUKernel.setArg(11, bufferV10);
	calcNUKernel.setArg(12, bufferQX);
	calcNUKernel.setArg(13, bufferQY);
	calcNUKernel.setArg(14, bufferUA);
	calcNUKernel.setArg(15, bufferVA);
	calcNUKernel.setArg(16, bufferNU);

	interpolateQUV10Kernel = cl::Kernel(program, "utils_interpolate_plane");

	interpolateQUV10Kernel.setArg(0, (int)data.ny);

	calcMaxPlaneKernel = cl::Kernel(program, "utils_max_plane");

	calcMaxPlaneKernel.setArg(0, (int)data.ny);
	calcMaxPlaneKernel.setArg(1, (int)data.nz + 1);
	calcMaxPlaneKernel.setArg(2, bufferNU);
	calcMaxPlaneKernel.setArg(3, bufferP);

	calcMaxRowKernel = cl::Kernel(program, "utils_max_row");

	calcMaxRowKernel.setArg(0, (int)data.ny);
	calcMaxRowKernel.setArg(1, bufferP);
	calcMaxRowKernel.setArg(2, bufferR);

	calcMaxKernel = cl::Kernel(program, "utils_max");

	calcMaxKernel.setArg(0, (int)data.nx);
	calcMaxKernel.setArg(1, bufferR);
	calcMaxKernel.setArg(2, bufferV);
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
	vector<string>& sourcePaths
) {
	for (auto& path : sourcePaths) {
		sources.push_back(
			loadKernelSource(
				path
			)
		);
	}
}

void Kernel::printDeviceDescription(cl::Platform& platform, cl::Device& device) {
	cout << "Platform: " << platform.getInfo<CL_PLATFORM_VERSION>() << endl;
	cout << "    Device: " << device.getInfo<CL_DEVICE_NAME>() << endl;
	cout << "    Vendor: " << device.getInfo<CL_DEVICE_VENDOR>() << endl;
	cout << "    Compute Units: " << device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << endl;
	cout << "    Clock Frequency: " << device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << " MHz" << endl;
	cout << "    Global Memory: " << device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() / (1024 * 1024) << " MB" << endl;
	cout << "    Local Memory: " << device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() / 1024 << " KB" << endl;
	cout << "    Max Work Group Size: " << device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << endl;
}

void Kernel::calcUA(double dt) {
	calcUAKernel.setArg(5, dt);
	calcUAKernel.setArg(12, bufferUA);
	calcUAKernel.setArg(13, bufferVA);
	calcUAKernel.setArg(17, bufferU1A);

	auto err = queue.enqueueNDRangeKernel(calcUAKernel, cl::NullRange, surfaceRange, cl::NullRange);	
}

void Kernel::calcVA(double dt) {
	calcVAKernel.setArg(5, dt);
	calcVAKernel.setArg(12, bufferUA);
	calcVAKernel.setArg(13, bufferVA);
	calcVAKernel.setArg(17, bufferV1A);

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

void Kernel::createTridiagonalMatrix(double dt) {
	createTridiagonalMatrixKernel.setArg(1, dt);

	auto err = queue.enqueueNDRangeKernel(createTridiagonalMatrixKernel, cl::NullRange, volumeRange, cl::NullRange);
}

void Kernel::calcResidual() {
	auto err = queue.enqueueNDRangeKernel(calcResidualKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

void Kernel::calcRHSU(double dt) {
	calcRHSUKernel.setArg(3, dt);
	calcRHSUKernel.setArg(10, bufferUA);

	auto err = queue.enqueueNDRangeKernel(calcRHSUKernel, cl::NullRange, volumeRange, cl::NullRange);

	err = queue.enqueueCopyBuffer(bufferUD, bufferUDT, 0, 0, sizeof(double) * data.nx * data.ny * data.nz);
}

void Kernel::calcUF() {
	auto err = queue.enqueueNDRangeKernel(calcUFKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

void Kernel::calcRHSV(double dt) {
	calcRHSVKernel.setArg(3, dt);
	calcRHSVKernel.setArg(10, bufferVA);

	auto err = queue.enqueueNDRangeKernel(calcRHSVKernel, cl::NullRange, volumeRange, cl::NullRange);

	err = queue.enqueueCopyBuffer(bufferVD, bufferVDT, 0, 0, sizeof(double) * data.nx * data.ny * data.nz);
}

void Kernel::calcVF() {
	auto err = queue.enqueueNDRangeKernel(calcVFKernel, cl::NullRange, surfaceRange, cl::NullRange);
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

void Kernel::updateQUV10Device(bool interpolateWind) {
	auto err = queue.enqueueWriteBuffer(bufferU10, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex].speed.u10.data());
	err = queue.enqueueWriteBuffer(bufferV10, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex].speed.v10.data());

	err = queue.enqueueWriteBuffer(bufferQX, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex].stress.qx.data());
	err = queue.enqueueWriteBuffer(bufferQY, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex].stress.qy.data());

	if (interpolateWind && data.currentWindIndex + 1 < data.windData.size()) {
		err = queue.enqueueCopyBuffer(bufferU10, bufferU10C, 0, 0, sizeof(double) * data.nx * data.ny);
		err = queue.enqueueCopyBuffer(bufferV10, bufferV10C, 0, 0, sizeof(double) * data.nx * data.ny);

		err = queue.enqueueCopyBuffer(bufferQX, bufferQXC, 0, 0, sizeof(double) * data.nx * data.ny);
		err = queue.enqueueCopyBuffer(bufferQY, bufferQYC, 0, 0, sizeof(double) * data.nx * data.ny);
		
		err = queue.enqueueWriteBuffer(bufferU10N, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex + 1].speed.u10.data());
		err = queue.enqueueWriteBuffer(bufferV10N, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex + 1].speed.v10.data());

		err = queue.enqueueWriteBuffer(bufferQXN, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex + 1].stress.qx.data());
		err = queue.enqueueWriteBuffer(bufferQYN, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, data.windData[data.currentWindIndex + 1].stress.qy.data());
	}	
}

void Kernel::interpolateQUV10(double alpha) {
	interpolateQUV10Kernel.setArg(1, alpha);	
	interpolateQUV10Kernel.setArg(2, bufferU10C);
	interpolateQUV10Kernel.setArg(3, bufferU10N);
	interpolateQUV10Kernel.setArg(4, bufferU10);

	auto err = queue.enqueueNDRangeKernel(interpolateQUV10Kernel, cl::NullRange, surfaceRange, cl::NullRange);
	
	interpolateQUV10Kernel.setArg(2, bufferV10C);
	interpolateQUV10Kernel.setArg(3, bufferV10N);
	interpolateQUV10Kernel.setArg(4, bufferV10);

	err = queue.enqueueNDRangeKernel(interpolateQUV10Kernel, cl::NullRange, surfaceRange, cl::NullRange);
	
	interpolateQUV10Kernel.setArg(2, bufferQXC);
	interpolateQUV10Kernel.setArg(3, bufferQXN);
	interpolateQUV10Kernel.setArg(4, bufferQX);

	err = queue.enqueueNDRangeKernel(interpolateQUV10Kernel, cl::NullRange, surfaceRange, cl::NullRange);
	
	interpolateQUV10Kernel.setArg(2, bufferQYC);
	interpolateQUV10Kernel.setArg(3, bufferQYN);
	interpolateQUV10Kernel.setArg(4, bufferQY);

	err = queue.enqueueNDRangeKernel(interpolateQUV10Kernel, cl::NullRange, surfaceRange, cl::NullRange);
}

void Kernel::updateQUV10(
	vector<double>& u10, vector<double>& v10,
	vector<double>& qx, vector<double>& qy
) {
	auto err = queue.enqueueReadBuffer(bufferU10, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, u10.data());
	err = queue.enqueueReadBuffer(bufferV10, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, v10.data());
	err = queue.enqueueReadBuffer(bufferQX, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, qx.data());
	err = queue.enqueueReadBuffer(bufferQY, CL_TRUE, 0, sizeof(double) * data.nx * data.ny, qy.data());
}

void Kernel::calcNU() {
	calcNUKernel.setArg(14, bufferUA);
	calcNUKernel.setArg(15, bufferVA);

	auto err = queue.enqueueNDRangeKernel(calcNUKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

double Kernel::calcNUMax() {
	double nuMax = 0;

	auto err = queue.enqueueNDRangeKernel(calcMaxPlaneKernel, cl::NullRange, surfaceRange, cl::NullRange);
	err = queue.enqueueNDRangeKernel(calcMaxRowKernel, cl::NullRange, rowRange, cl::NullRange);
	err = queue.enqueueNDRangeKernel(calcMaxKernel, cl::NullRange, valueRange, cl::NullRange);

	err = queue.enqueueReadBuffer(bufferV, CL_TRUE, 0, sizeof(double), &nuMax);

	return nuMax;
}

void Kernel::updateNU() {
	auto err = queue.enqueueReadBuffer(bufferNU, CL_TRUE, 0, sizeof(double) * data.nx * data.ny * (data.nz + 1), data.nu.data());
}

path Kernel::addDirectory(path outDir) {
	return outDir / path("dp");
}

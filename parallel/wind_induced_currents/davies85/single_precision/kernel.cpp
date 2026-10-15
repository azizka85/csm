#include <iostream>
#include <fstream>

#include <stdexcept>

#include "kernel.h"

using namespace Parallel::WindInducedCurrents::Davies85::SinglePrecision;

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

	vector<float> dzF(data.dz.begin(), data.dz.end());
	bufferDZ = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nz);
	err = queue.enqueueWriteBuffer(bufferDZ, CL_TRUE, 0, sizeof(float) * data.nz, dzF.data());

	vector<float> hF(data.bathymetry.h.begin(), data.bathymetry.h.end());
	bufferH = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferH, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, hF.data());

	vector<float> uaF(data.ua.begin(), data.ua.end());
	bufferUA = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferUA, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, uaF.data());

	bufferU1A = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);

	vector<float> vaF(data.va.begin(), data.va.end());
	bufferVA = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferVA, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, vaF.data());

	bufferV1A = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);

	vector<float> zF(data.z.begin(), data.z.end());
	bufferZ = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferZ, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, zF.data());

	vector<float> ufF(data.uf.begin(), data.uf.end());
	bufferUF = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny * data.nz);
	err = queue.enqueueWriteBuffer(bufferUF, CL_TRUE, 0, sizeof(float) * data.nx * data.ny * data.nz, ufF.data());

	bufferUD = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny * data.nz);
	bufferUDT = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny * data.nz);

	vector<float> vfF(data.vf.begin(), data.vf.end());
	bufferVF = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny * data.nz);
	err = queue.enqueueWriteBuffer(bufferVF, CL_TRUE, 0, sizeof(float) * data.nx * data.ny * data.nz, vfF.data());

	bufferVD = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny * data.nz);
	bufferVDT = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny * data.nz);

	bufferAL = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny * data.nz);
	bufferAC = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny * data.nz);
	bufferAR = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny * data.nz);

	vector<float> u10F(data.windData[data.currentWindIndex].speed.u10.begin(), data.windData[data.currentWindIndex].speed.u10.end());
	bufferU10 = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferU10, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, u10F.data());

	vector<float> v10F(data.windData[data.currentWindIndex].speed.v10.begin(), data.windData[data.currentWindIndex].speed.v10.end());
	bufferV10 = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferV10, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, v10F.data());

	vector<float> qxF(data.windData[data.currentWindIndex].stress.qx.begin(), data.windData[data.currentWindIndex].stress.qx.end());
	bufferQX = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferQX, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, qxF.data());

	vector<float> qyF(data.windData[data.currentWindIndex].stress.qy.begin(), data.windData[data.currentWindIndex].stress.qy.end());
	bufferQY = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
	err = queue.enqueueWriteBuffer(bufferQY, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, qyF.data());

	if (interpolateWind) {
		bufferU10C = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
		err = queue.enqueueCopyBuffer(bufferU10, bufferU10C, 0, 0, sizeof(float) * data.nx * data.ny);

		bufferV10C = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
		err = queue.enqueueCopyBuffer(bufferV10, bufferV10C, 0, 0, sizeof(float) * data.nx * data.ny);

		bufferQXC = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
		err = queue.enqueueCopyBuffer(bufferQX, bufferQXC, 0, 0, sizeof(float) * data.nx * data.ny);

		bufferQYC = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
		err = queue.enqueueCopyBuffer(bufferQY, bufferQYC, 0, 0, sizeof(float) * data.nx * data.ny);

		bufferU10N = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
		bufferV10N = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
		bufferQXN = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);
		bufferQYN = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);

		if (data.currentWindIndex + 1 < data.windData.size()) {
			vector<float> u10NF(data.windData[data.currentWindIndex + 1].speed.u10.begin(), data.windData[data.currentWindIndex + 1].speed.u10.end());
			err = queue.enqueueWriteBuffer(bufferU10N, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, u10NF.data());

			vector<float> v10NF(data.windData[data.currentWindIndex + 1].speed.v10.begin(), data.windData[data.currentWindIndex + 1].speed.v10.end());
			err = queue.enqueueWriteBuffer(bufferV10N, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, v10NF.data());

			vector<float> qxNF(data.windData[data.currentWindIndex + 1].stress.qx.begin(), data.windData[data.currentWindIndex + 1].stress.qx.end());
			err = queue.enqueueWriteBuffer(bufferQXN, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, qxNF.data());

			vector<float> qyNF(data.windData[data.currentWindIndex + 1].stress.qy.begin(), data.windData[data.currentWindIndex + 1].stress.qy.end());
			err = queue.enqueueWriteBuffer(bufferQYN, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, qyNF.data());
		}
		else {
			err = queue.enqueueCopyBuffer(bufferU10, bufferU10N, 0, 0, sizeof(float) * data.nx * data.ny);
			err = queue.enqueueCopyBuffer(bufferV10, bufferV10N, 0, 0, sizeof(float) * data.nx * data.ny);
			err = queue.enqueueCopyBuffer(bufferQX, bufferQXN, 0, 0, sizeof(float) * data.nx * data.ny);
			err = queue.enqueueCopyBuffer(bufferQY, bufferQYN, 0, 0, sizeof(float) * data.nx * data.ny);
		}
	}

	vector<float> nuF(data.nu.begin(), data.nu.end());
	bufferNU = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny * (data.nz + 1));
	err = queue.enqueueWriteBuffer(bufferNU, CL_TRUE, 0, sizeof(float) * data.nx * data.ny * (data.nz + 1), nuF.data());

	bufferV = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float));
	bufferR = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx);
	bufferP = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * data.nx * data.ny);

	calcUAKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_ua");

	calcUAKernel.setArg(0, (float)data.f);
	calcUAKernel.setArg(1, (float)data.kb);
	calcUAKernel.setArg(2, (float)data.g);
	calcUAKernel.setArg(3, (float)data.rho);
	calcUAKernel.setArg(4, (float)data.dx);
	calcUAKernel.setArg(6, (int)data.nx);
	calcUAKernel.setArg(7, (int)data.ny);
	calcUAKernel.setArg(8, (int)data.nz);
	calcUAKernel.setArg(9, bufferH);
	calcUAKernel.setArg(10, bufferQX);
	calcUAKernel.setArg(11, bufferUA);
	calcUAKernel.setArg(12, bufferVA);
	calcUAKernel.setArg(13, bufferZ);
	calcUAKernel.setArg(14, bufferUF);
	calcUAKernel.setArg(15, bufferU1A);

	calcVAKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_va");

	calcVAKernel.setArg(0, (float)data.f);
	calcVAKernel.setArg(1, (float)data.kb);
	calcVAKernel.setArg(2, (float)data.g);
	calcVAKernel.setArg(3, (float)data.rho);
	calcVAKernel.setArg(4, (float)data.dy);
	calcVAKernel.setArg(6, (int)data.nx);
	calcVAKernel.setArg(7, (int)data.ny);
	calcVAKernel.setArg(8, (int)data.nz);
	calcVAKernel.setArg(9, bufferH);
	calcVAKernel.setArg(10, bufferQY);
	calcVAKernel.setArg(11, bufferUA);
	calcVAKernel.setArg(12, bufferVA);
	calcVAKernel.setArg(13, bufferZ);
	calcVAKernel.setArg(14, bufferVF);
	calcVAKernel.setArg(15, bufferV1A);

	calcRHSKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_rhs");

	calcRHSKernel.setArg(0, (float)data.f);
	calcRHSKernel.setArg(1, (float)data.kb);
	calcRHSKernel.setArg(2, (float)data.rho);
	calcRHSKernel.setArg(3, (float)data.dx);
	calcRHSKernel.setArg(4, (float)data.dy);
	calcRHSKernel.setArg(6, (int)data.ny);
	calcRHSKernel.setArg(7, (int)data.nz);
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

	createTridiagonalMatrixKernel.setArg(0, (float)data.kb);
	createTridiagonalMatrixKernel.setArg(2, (int)data.ny);
	createTridiagonalMatrixKernel.setArg(3, (int)data.nz);
	createTridiagonalMatrixKernel.setArg(4, bufferNU);
	createTridiagonalMatrixKernel.setArg(5, bufferDZ);
	createTridiagonalMatrixKernel.setArg(6, bufferH);
	createTridiagonalMatrixKernel.setArg(7, bufferAL);
	createTridiagonalMatrixKernel.setArg(8, bufferAC);
	createTridiagonalMatrixKernel.setArg(9, bufferAR);

	calcUVFKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_uvf");

	calcUVFKernel.setArg(0, (int)data.ny);
	calcUVFKernel.setArg(1, (int)data.nz);
	calcUVFKernel.setArg(2, bufferH);
	calcUVFKernel.setArg(3, bufferAL);
	calcUVFKernel.setArg(4, bufferAC);
	calcUVFKernel.setArg(5, bufferAR);
	calcUVFKernel.setArg(6, bufferUF);
	calcUVFKernel.setArg(7, bufferVF);
	calcUVFKernel.setArg(8, bufferUD);
	calcUVFKernel.setArg(9, bufferVD);

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

	calcZKernel.setArg(0, (float)data.dx);
	calcZKernel.setArg(1, (float)data.dy);
	calcZKernel.setArg(3, (int)data.nx);
	calcZKernel.setArg(4, (int)data.ny);
	calcZKernel.setArg(5, (int)data.nz);
	calcZKernel.setArg(6, bufferH);
	calcZKernel.setArg(7, bufferUA);
	calcZKernel.setArg(8, bufferVA);
	calcZKernel.setArg(9, bufferZ);

	calcNUKernel = cl::Kernel(program, "viscosity");

	calcNUKernel.setArg(0, (float)data.k0);
	calcNUKernel.setArg(1, (float)data.rho);
	calcNUKernel.setArg(2, (float)data.f);
	calcNUKernel.setArg(3, (float)data.k);
	calcNUKernel.setArg(4, (float)data.sigma);
	calcNUKernel.setArg(5, (float)data.num);
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

void Kernel::loadKernelSources(cl::Program::Sources& sources, vector<string>& sourcePaths) {
	for (auto& path : sourcePaths) {
		sources.push_back(
			loadKernelSource(
				path
			)
		);
	}
}

void Kernel::readFloatBuffer(const cl::Buffer& buf, vector<double>& dst, size_t count) {
	if (scratch.size() < count) scratch.resize(count);

	queue.enqueueReadBuffer(buf, CL_TRUE, 0, sizeof(float) * count, scratch.data());
	copy(scratch.begin(), scratch.begin() + count, dst.begin());
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
	calcUAKernel.setArg(5, (float)dt);
	calcUAKernel.setArg(11, bufferUA);
	calcUAKernel.setArg(12, bufferVA);
	calcUAKernel.setArg(15, bufferU1A);

	auto err = queue.enqueueNDRangeKernel(calcUAKernel, cl::NullRange, surfaceRange, cl::NullRange);	
}

void Kernel::calcVA(double dt) {
	calcVAKernel.setArg(5, (float)dt);
	calcVAKernel.setArg(11, bufferUA);
	calcVAKernel.setArg(12, bufferVA);
	calcVAKernel.setArg(15, bufferV1A);

	auto err = queue.enqueueNDRangeKernel(calcVAKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

void Kernel::updateUVA() {
	readFloatBuffer(bufferUA, data.ua, data.nx * data.ny);
	readFloatBuffer(bufferVA, data.va, data.nx * data.ny);
}

void Kernel::swapUA() {
	swap(bufferU1A, bufferUA);
}

void Kernel::swapVA() {
	swap(bufferV1A, bufferVA);
}

void Kernel::calcRHS(double dt) {
	calcRHSKernel.setArg(5, (float)dt);
	calcRHSKernel.setArg(13, bufferUA);
	calcRHSKernel.setArg(14, bufferVA);

	auto err = queue.enqueueNDRangeKernel(calcRHSKernel, cl::NullRange, volumeRange, cl::NullRange);

	err = queue.enqueueCopyBuffer(bufferUD, bufferUDT, 0, 0, sizeof(float) * data.nx * data.ny * data.nz);
	err = queue.enqueueCopyBuffer(bufferVD, bufferVDT, 0, 0, sizeof(float) * data.nx * data.ny * data.nz);
}

void Kernel::createTridiagonalMatrix(double dt) {
	createTridiagonalMatrixKernel.setArg(1, (float)dt);

	auto err = queue.enqueueNDRangeKernel(createTridiagonalMatrixKernel, cl::NullRange, volumeRange, cl::NullRange);
}

void Kernel::calcResidual() {
	auto err = queue.enqueueNDRangeKernel(calcResidualKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

void Kernel::calcUVF() {
	auto err = queue.enqueueNDRangeKernel(calcUVFKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

void Kernel::updateUVF() {
	readFloatBuffer(bufferUF, data.uf, data.nx * data.ny * data.nz);
	readFloatBuffer(bufferVF, data.vf, data.nx * data.ny * data.nz);
}

void Kernel::updateUVD() {
	readFloatBuffer(bufferUDT, data.ud, data.nx * data.ny * data.nz);
	readFloatBuffer(bufferVDT, data.vd, data.nx * data.ny * data.nz);
}

void Kernel::calcZ(double dt) {
	calcZKernel.setArg(2, (float)dt);
	calcZKernel.setArg(7, bufferUA);
	calcZKernel.setArg(8, bufferVA);

	auto err = queue.enqueueNDRangeKernel(calcZKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

void Kernel::updateZ() {
	readFloatBuffer(bufferZ, data.z, data.nx * data.ny);
}

void Kernel::updateQUV10Device(bool interpolateWind) {
	vector<float> u10F(data.windData[data.currentWindIndex].speed.u10.begin(), data.windData[data.currentWindIndex].speed.u10.end());
	auto err = queue.enqueueWriteBuffer(bufferU10, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, u10F.data());

	vector<float> v10F(data.windData[data.currentWindIndex].speed.v10.begin(), data.windData[data.currentWindIndex].speed.v10.end());
	err = queue.enqueueWriteBuffer(bufferV10, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, v10F.data());

	vector<float> qxF(data.windData[data.currentWindIndex].stress.qx.begin(), data.windData[data.currentWindIndex].stress.qx.end());
	err = queue.enqueueWriteBuffer(bufferQX, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, qxF.data());

	vector<float> qyF(data.windData[data.currentWindIndex].stress.qy.begin(), data.windData[data.currentWindIndex].stress.qy.end());
	err = queue.enqueueWriteBuffer(bufferQY, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, qyF.data());

	if (interpolateWind && data.currentWindIndex + 1 < data.windData.size()) {
		err = queue.enqueueCopyBuffer(bufferU10, bufferU10C, 0, 0, sizeof(double) * data.nx * data.ny);
		err = queue.enqueueCopyBuffer(bufferV10, bufferV10C, 0, 0, sizeof(double) * data.nx * data.ny);

		err = queue.enqueueCopyBuffer(bufferQX, bufferQXC, 0, 0, sizeof(double) * data.nx * data.ny);
		err = queue.enqueueCopyBuffer(bufferQY, bufferQYC, 0, 0, sizeof(double) * data.nx * data.ny);

		vector<float> u10NF(data.windData[data.currentWindIndex + 1].speed.u10.begin(), data.windData[data.currentWindIndex + 1].speed.u10.end());
		err = queue.enqueueWriteBuffer(bufferU10N, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, u10NF.data());

		vector<float> v10NF(data.windData[data.currentWindIndex + 1].speed.v10.begin(), data.windData[data.currentWindIndex + 1].speed.v10.end());
		err = queue.enqueueWriteBuffer(bufferV10N, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, v10NF.data());

		vector<float> qxNF(data.windData[data.currentWindIndex + 1].stress.qx.begin(), data.windData[data.currentWindIndex + 1].stress.qx.end());
		err = queue.enqueueWriteBuffer(bufferQXN, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, qxNF.data());

		vector<float> qyNF(data.windData[data.currentWindIndex + 1].stress.qy.begin(), data.windData[data.currentWindIndex + 1].stress.qy.end());
		err = queue.enqueueWriteBuffer(bufferQYN, CL_TRUE, 0, sizeof(float) * data.nx * data.ny, qyNF.data());
	}
}

void Kernel::interpolateQUV10(double alpha) {
	interpolateQUV10Kernel.setArg(1, (float)alpha);
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
	readFloatBuffer(bufferU10, u10, data.nx * data.ny);
	readFloatBuffer(bufferV10, v10, data.nx * data.ny);
	readFloatBuffer(bufferQX, qx, data.nx * data.ny);
	readFloatBuffer(bufferQY, qy, data.nx * data.ny);
}

void Kernel::calcNU() {
	calcNUKernel.setArg(14, bufferUA);
	calcNUKernel.setArg(15, bufferVA);

	auto err = queue.enqueueNDRangeKernel(calcNUKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

double Kernel::calcNUMax() {
	float nuMax = 0;

	auto err = queue.enqueueNDRangeKernel(calcMaxPlaneKernel, cl::NullRange, surfaceRange, cl::NullRange);
	err = queue.enqueueNDRangeKernel(calcMaxRowKernel, cl::NullRange, rowRange, cl::NullRange);
	err = queue.enqueueNDRangeKernel(calcMaxKernel, cl::NullRange, valueRange, cl::NullRange);

	err = queue.enqueueReadBuffer(bufferV, CL_TRUE, 0, sizeof(float), &nuMax);

	return nuMax;
}

void Kernel::updateNU() {
	readFloatBuffer(bufferNU, data.nu, data.nx * data.ny * (data.nz + 1));
}

path Kernel::addDirectory(path outDir) {
	return outDir / path("sp");
}

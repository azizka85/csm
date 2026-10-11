#include <fstream>

#include <stdexcept>

#include "kernel.h"

using namespace Parallel::SinglePrecision::WindInducedCurrents::Davies85;

Kernel::Kernel(
	int nx, int ny,
	double dx, double dy,
	double rho, double f, double kb, double g,
	vector<double>& dz, vector<double>& h,
	vector<double>& ua, vector<double>& va,
	vector<double>& z,
	vector<double>& uf, vector<double>& vf,
	vector<double>& ud, vector<double>& vd,
	vector<double>& u10, vector<double>& v10,
	vector<double>& qx, vector<double>& qy,
	vector<double>& nu
): nx(nx), ny(ny),
	dx(dx), dy(dy),
	rho(rho), f(f), kb(kb), g(g),
	dz(dz), h(h),
	ua(ua), va(va),
	z(z),
	uf(uf), vf(vf),
	ud(ud), vd(vd),
	u10(u10), v10(v10),
	qx(qx), qy(qy),
	nu(nu) 
{
	nz = dz.size();

	surfaceRange = cl::NDRange(nx, ny);
	volumeRange = cl::NDRange(nx, ny, nz);

	vector<cl::Platform> platforms;

	cl::Platform::get(&platforms);

	cl::Platform platform = platforms.front();

	vector<cl::Device> devices;

	platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);

	cl::Device device = devices.front();
	cl::Context context(device);
	queue = cl::CommandQueue(context, device);

	cl::Program::Sources sources;

	loadKernelSources(sources);

	cl::Program program(context, sources);

	auto err = program.build(device);

	if (err == -11) {
		auto build_log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>();

		bool ok = true;
	}

	vector<float> dzF(dz.begin(), dz.end());
	bufferDZ = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nz);
	err = queue.enqueueWriteBuffer(bufferDZ, CL_TRUE, 0, sizeof(float) * nz, dzF.data());

	vector<float> hF(h.begin(), h.end());
	bufferH = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
	err = queue.enqueueWriteBuffer(bufferH, CL_TRUE, 0, sizeof(float) * nx * ny, hF.data());

	vector<float> uaF(ua.begin(), ua.end());
	bufferUA = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
	err = queue.enqueueWriteBuffer(bufferUA, CL_TRUE, 0, sizeof(float) * nx * ny, uaF.data());

	bufferU1A = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);

	vector<float> vaF(va.begin(), va.end());
	bufferVA = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
	err = queue.enqueueWriteBuffer(bufferVA, CL_TRUE, 0, sizeof(float) * nx * ny, vaF.data());

	bufferV1A = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);

	vector<float> zF(z.begin(), z.end());
	bufferZ = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
	err = queue.enqueueWriteBuffer(bufferZ, CL_TRUE, 0, sizeof(float) * nx * ny, zF.data());

	vector<float> ufF(uf.begin(), uf.end());
	bufferUF = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);
	err = queue.enqueueWriteBuffer(bufferUF, CL_TRUE, 0, sizeof(float) * nx * ny * nz, ufF.data());

	bufferUD = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);
	bufferUDT = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);

	vector<float> vfF(vf.begin(), vf.end());
	bufferVF = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);
	err = queue.enqueueWriteBuffer(bufferVF, CL_TRUE, 0, sizeof(float) * nx * ny * nz, vfF.data());

	bufferVD = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);
	bufferVDT = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);

	bufferAL = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);
	bufferAC = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);
	bufferAR = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);

	vector<float> u10F(u10.begin(), u10.end());
	bufferU10 = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
	err = queue.enqueueWriteBuffer(bufferU10, CL_TRUE, 0, sizeof(float) * nx * ny, u10F.data());

	vector<float> v10F(v10.begin(), v10.end());
	bufferV10 = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
	err = queue.enqueueWriteBuffer(bufferV10, CL_TRUE, 0, sizeof(float) * nx * ny, v10F.data());

	vector<float> qxF(qx.begin(), qx.end());
	bufferQX = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
	err = queue.enqueueWriteBuffer(bufferQX, CL_TRUE, 0, sizeof(float) * nx * ny, qxF.data());

	vector<float> qyF(qy.begin(), qy.end());
	bufferQY = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
	err = queue.enqueueWriteBuffer(bufferQY, CL_TRUE, 0, sizeof(float) * nx * ny, qyF.data());

	vector<float> nuF(nu.begin(), nu.end());
	bufferNU = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * (nz + 1));
	err = queue.enqueueWriteBuffer(bufferNU, CL_TRUE, 0, sizeof(float) * nx * ny * (nz + 1), nuF.data());

	calcUAKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_ua");

	calcUAKernel.setArg(0, (float)f);
	calcUAKernel.setArg(1, (float)kb);
	calcUAKernel.setArg(2, (float)g);
	calcUAKernel.setArg(3, (float)rho);
	calcUAKernel.setArg(4, (float)dx);
	calcUAKernel.setArg(6, nx);
	calcUAKernel.setArg(7, ny);
	calcUAKernel.setArg(8, nz);
	calcUAKernel.setArg(9, bufferH);
	calcUAKernel.setArg(10, bufferQX);
	calcUAKernel.setArg(11, bufferUA);
	calcUAKernel.setArg(12, bufferVA);
	calcUAKernel.setArg(13, bufferZ);
	calcUAKernel.setArg(14, bufferUF);
	calcUAKernel.setArg(15, bufferU1A);

	calcVAKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_va");

	calcVAKernel.setArg(0, (float)f);
	calcVAKernel.setArg(1, (float)kb);
	calcVAKernel.setArg(2, (float)g);
	calcVAKernel.setArg(3, (float)rho);
	calcVAKernel.setArg(4, (float)dy);
	calcVAKernel.setArg(6, nx);
	calcVAKernel.setArg(7, ny);
	calcVAKernel.setArg(8, nz);
	calcVAKernel.setArg(9, bufferH);
	calcVAKernel.setArg(10, bufferQY);
	calcVAKernel.setArg(11, bufferUA);
	calcVAKernel.setArg(12, bufferVA);
	calcVAKernel.setArg(13, bufferZ);
	calcVAKernel.setArg(14, bufferVF);
	calcVAKernel.setArg(15, bufferV1A);

	calcRHSKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_rhs");

	calcRHSKernel.setArg(0, (float)f);
	calcRHSKernel.setArg(1, (float)kb);
	calcRHSKernel.setArg(2, (float)rho);
	calcRHSKernel.setArg(3, (float)dx);
	calcRHSKernel.setArg(4, (float)dy);
	calcRHSKernel.setArg(6, ny);
	calcRHSKernel.setArg(7, nz);
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

	createTridiagonalMatrixKernel.setArg(0, (float)kb);
	createTridiagonalMatrixKernel.setArg(2, ny);
	createTridiagonalMatrixKernel.setArg(3, nz);
	createTridiagonalMatrixKernel.setArg(4, bufferNU);
	createTridiagonalMatrixKernel.setArg(5, bufferDZ);
	createTridiagonalMatrixKernel.setArg(6, bufferH);
	createTridiagonalMatrixKernel.setArg(7, bufferAL);
	createTridiagonalMatrixKernel.setArg(8, bufferAC);
	createTridiagonalMatrixKernel.setArg(9, bufferAR);

	calcUVFKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_uvf");

	calcUVFKernel.setArg(0, ny);
	calcUVFKernel.setArg(1, nz);
	calcUVFKernel.setArg(2, bufferH);
	calcUVFKernel.setArg(3, bufferAL);
	calcUVFKernel.setArg(4, bufferAC);
	calcUVFKernel.setArg(5, bufferAR);
	calcUVFKernel.setArg(6, bufferUF);
	calcUVFKernel.setArg(7, bufferVF);
	calcUVFKernel.setArg(8, bufferUD);
	calcUVFKernel.setArg(9, bufferVD);

	calcResidualKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_residual");

	calcResidualKernel.setArg(0, ny);
	calcResidualKernel.setArg(1, nz);
	calcResidualKernel.setArg(2, bufferH);
	calcResidualKernel.setArg(3, bufferAL);
	calcResidualKernel.setArg(4, bufferAC);
	calcResidualKernel.setArg(5, bufferAR);
	calcResidualKernel.setArg(6, bufferUF);
	calcResidualKernel.setArg(7, bufferVF);
	calcResidualKernel.setArg(8, bufferUDT);
	calcResidualKernel.setArg(9, bufferVDT);

	calcZKernel = cl::Kernel(program, "wind_induced_currents_davies85_calc_z");

	calcZKernel.setArg(0, (float)dx);
	calcZKernel.setArg(1, (float)dy);
	calcZKernel.setArg(3, nx);
	calcZKernel.setArg(4, ny);
	calcZKernel.setArg(5, nz);
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

void Kernel::loadKernelSources(cl::Program::Sources& sources) {
	sources.push_back(
		loadKernelSource(
			"opencl/single_precision/wind_induced_currents/davies85/solver.cl"
		)
	);
}

void Kernel::readFloatBuffer(const cl::Buffer& buf, vector<double>& dst, size_t count) {
	if (scratch.size() < count) scratch.resize(count);

	queue.enqueueReadBuffer(buf, CL_TRUE, 0, sizeof(float) * count, scratch.data());
	copy(scratch.begin(), scratch.begin() + count, dst.begin());
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
	readFloatBuffer(bufferUA, ua, nx * ny);
	readFloatBuffer(bufferVA, va, nx * ny);
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

	err = queue.enqueueCopyBuffer(bufferUD, bufferUDT, 0, 0, sizeof(float) * nx * ny * nz);
	err = queue.enqueueCopyBuffer(bufferVD, bufferVDT, 0, 0, sizeof(float) * nx * ny * nz);
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
	readFloatBuffer(bufferUF, uf, nx * ny * nz);
	readFloatBuffer(bufferVF, vf, nx * ny * nz);
}

void Kernel::updateUVD() {
	readFloatBuffer(bufferUDT, ud, nx * ny * nz);
	readFloatBuffer(bufferVDT, vd, nx * ny * nz);
}

void Kernel::calcZ(double dt) {
	calcZKernel.setArg(2, (float)dt);
	calcZKernel.setArg(7, bufferUA);
	calcZKernel.setArg(8, bufferVA);

	auto err = queue.enqueueNDRangeKernel(calcZKernel, cl::NullRange, surfaceRange, cl::NullRange);
}

void Kernel::updateZ() {
	readFloatBuffer(bufferZ, z, nx * ny);
}

void Kernel::updateQUV10Device() {
	vector<float> u10F(u10.begin(), u10.end());
	auto err = queue.enqueueWriteBuffer(bufferU10, CL_TRUE, 0, sizeof(float) * nx * ny, u10F.data());

	vector<float> v10F(v10.begin(), v10.end());
	err = queue.enqueueWriteBuffer(bufferV10, CL_TRUE, 0, sizeof(float) * nx * ny, v10F.data());

	vector<float> qxF(qx.begin(), qx.end());
	err = queue.enqueueWriteBuffer(bufferQX, CL_TRUE, 0, sizeof(float) * nx * ny, qxF.data());

	vector<float> qyF(qy.begin(), qy.end());
	err = queue.enqueueWriteBuffer(bufferQY, CL_TRUE, 0, sizeof(float) * nx * ny, qyF.data());
}

void Kernel::calcNU() {}

void Kernel::updateNU() {}




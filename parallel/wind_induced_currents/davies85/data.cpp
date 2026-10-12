#include "data.h"

using namespace Parallel::WindInducedCurrents;

Davies85::Data::Data(
	double g, double kb,
	Calc::Data::PhysicalParams phys,	
	double dx, double dy,
	shared_ptr<Generators::Area::IGenerator> areaGenerator,
	shared_ptr<Generators::DZ::IGenerator> dzGenerator,
	shared_ptr<Generators::Bathymetry::IGenerator> hGenerator,
	shared_ptr<Generators::Wind::IGenerator> qGenerator
) {
	this->g = g;
	this->kb = kb;

	rho = phys.rho;
	f = phys.f;
	
	this->dx = dx;
	this->dy = dy;

	this->areaGenerator = areaGenerator;
	this->dzGenerator = dzGenerator;
	this->hGenerator = hGenerator;
	this->qGenerator = qGenerator;
}

void Davies85::Data::init() {
	auto geo = areaGenerator->generate();

	nx = static_cast<size_t>(
		ceil(geo.l / dx)
	) + 1;

	ny = static_cast<size_t>(
		ceil(geo.w / dy)
	) + 1;

	dz = dzGenerator->generate();

	nz = dz.size();

	h = hGenerator->generate(nx, ny);

	ua = vector<double>(nx * ny);
	va = vector<double>(nx * ny);

	z = vector<double>(nx * ny);

	uf = vector<double>(nx * ny * nz);
	vf = vector<double>(nx * ny * nz);

	ud = vector<double>(nx * ny * nz);
	vd = vector<double>(nx * ny * nz);

	nu = vector<double>(nx * ny * (nz + 1));

	currentWindIndex = 0;

	windData = qGenerator->generate(nx, ny);
}

path Davies85::Data::addDirectory(path outDir) {
	outDir = hGenerator->addDirectory(outDir);
	outDir = qGenerator->addDirectory(outDir);
	outDir = areaGenerator->addDirectory(outDir);

	outDir = outDir / format("dx={}, dy={}", dx, dy);

	outDir = dzGenerator->addDirectory(outDir);

	return outDir;
}

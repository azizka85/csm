#include <calc/data.h>

#include <utils/data.h>

#include "data.h"

using namespace Parallel::WindInducedCurrents;

Davies85::Data::Data(
	double g, double kb, double k0, double k, double sigma, double num,
	Calc::Data::PhysicalParams phys,	
	double dx, double dy,
	shared_ptr<Generators::Area::IGenerator> areaGenerator,
	shared_ptr<Generators::DZ::IGenerator> dzGenerator,
	shared_ptr<Generators::Bathymetry::IGenerator> hGenerator,
	shared_ptr<Generators::Wind::IGenerator> qGenerator,
	shared_ptr<Generators::Viscosity::IGenerator> nuGenerator
) {
	this->g = g;
	this->kb = kb;
	this->k0 = k0;
	this->k = k;
	this->sigma = sigma;

	this->num = num;

	rho = phys.rho;
	f = phys.f;
	
	this->dx = dx;
	this->dy = dy;

	this->areaGenerator = areaGenerator;
	this->dzGenerator = dzGenerator;
	this->hGenerator = hGenerator;
	this->qGenerator = qGenerator;
	this->nuGenerator = nuGenerator;
}

void Davies85::Data::init(bool symmetrise) {
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

	windData = qGenerator->generate(nx, ny);

	if (symmetrise) {
		h = Utils::Data::symmetrise(nx, ny, h);

		for (auto& wind : windData) {
			wind.speed.u10 = Utils::Data::symmetrise(nx, ny, wind.speed.u10);
			wind.speed.v10 = Utils::Data::symmetrise(nx, ny, wind.speed.v10);

			wind.stress.qx = Utils::Data::symmetrise(nx, ny, wind.stress.qx);
			wind.stress.qy = Utils::Data::symmetrise(nx, ny, wind.stress.qy);
		}

		nx *= 2;
		ny *= 2;
	}

	ua = vector<double>(nx * ny);
	va = vector<double>(nx * ny);

	z = vector<double>(nx * ny);

	uf = vector<double>(nx * ny * nz);
	vf = vector<double>(nx * ny * nz);

	ud = vector<double>(nx * ny * nz);
	vd = vector<double>(nx * ny * nz);

	currentWindIndex = 0;

	nu = nuGenerator->generate(
		Calc::Data::VolumeState{
			.nx = nx, .ny = ny, .nz = nz + 1,
			.dx = dx, .dy = dy,
			.dz = dz, .h = h
		},
		Calc::Data::WindCurrentData{
			.speed = Calc::Data::WindSpeedData {
				.u10 = windData[currentWindIndex].speed.u10,
				.v10 = windData[currentWindIndex].speed.v10
			},
			.stress = Calc::Data::WindStressData {
				.qx = windData[currentWindIndex].stress.qx,
				.qy = windData[currentWindIndex].stress.qy
			},
			.current = Calc::Data::Current {
				.ua = ua,
				.va = va
			}
		}
	);
}

path Davies85::Data::addDirectory(path outDir) {
	outDir = hGenerator->addDirectory(outDir);
	outDir = qGenerator->addDirectory(outDir);
	outDir = areaGenerator->addDirectory(outDir);

	outDir = outDir / format("dx={}, dy={}", dx, dy);

	outDir = dzGenerator->addDirectory(outDir);

	return outDir;
}

#include <utils/checks/fs.h>

#include "surface.h"

using namespace Utils::Checks;

using namespace Writers::Surface;

Writer::Writer(string fileName, WriteSurfaceDataHeaderFn writeSurfaceDataHeader) {
	this->fileName = fileName;
	this->writeSurfaceDataHeader = writeSurfaceDataHeader;
}

void Writer::writeSurfaceVector(string name, size_t precision, SurfaceState surface, SurfaceVector vector, ofstream file) {
	auto nx = surface.nx;
	auto ny = surface.ny;

	auto vx = vector.vx;
	auto vy = vector.vy;

	file << format("VECTORS {} double", name) << endl;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			int p = j + i * ny;

			file << format("{:.{}f} {:.{}f} 0", vx[p], precision, vy[p], precision) << endl;
		}
	}
}

void Writer::writeSurfaceScalar(string name, size_t precision, SurfaceState surface, vector<double> data, ofstream file) {
	auto nx = surface.nx;
	auto ny = surface.ny;

	file << format("SCALARS {} double", name) << endl;
	file << "LOOKUP_TABLE default" << endl;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			int p = j + i * ny;

			file << format("{:.{}f}", data[p], precision) << endl;
		}
	}
}

template<size_t N>
inline path Writer::write(
	TimeState time, 
	SurfaceState surface, 
	span<Column, N> columns, 
	span<DataKind, N> data, 
	path outDir
) {
	FS::checkPath(outDir);

	auto fName = fileName.empty() ? format("data.{:03}.vtk", time.m) : fileName;

	auto filePath = outDir / path(fileName);

	ofstream file(filePath);

	writeSurfaceDataHeader(time, surface, file);

	for (size_t i = 0; i < N; i++) {
		if (holds_alternative<SurfaceVector>(data[i])) {
			writeSurfaceVector(columns[i].name, columns[i].precision, surface, data[i], file);
		}
		else if (holds_alternative<vector<double>>(data[i])) {
			writeSurfaceScalar(columns[i].name, columns[i].precision, surface, data[i], file);
		}
	}

	return filePath;
}

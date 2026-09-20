#include <format>
#include <utils/checks/fs.h>

#include "volume.h"

using namespace Utils::Checks;

using namespace Writers::Volume;

Writer::Writer(string fileName, WriteVolumeDataHeaderFn writeVolumeDataHeader) {
	this->fileName = fileName;
	this->writeVolumeDataHeader = writeVolumeDataHeader;
}

path Writer::write(TimeState time, VolumeState volume, const vector<Column> &columns, const vector<DataKind> &data, path outDir) {
	FS::checkPath(outDir);

	auto fName = fileName.empty() ? format("data.{:03}.vtk", time.m) : fileName;

	auto filePath = outDir / path(fName);

	ofstream file(filePath);

	writeVolumeDataHeader(time, volume, file);

	for (size_t i = 0; i < columns.size() && i < data.size(); i++) {
		if (holds_alternative<VolumeVector>(data[i])) {
			writeVolumeVector(columns[i].name, columns[i].precision, volume, get<VolumeVector>(data[i]), file);
		}
		else if (holds_alternative<vector<double>>(data[i])) {
			writeVolumeScalar(columns[i].name, columns[i].precision, volume, get<vector<double>>(data[i]), file);
		}
	}

	return filePath;
}

void Writer::writeVolumeVector(string name, size_t precision, VolumeState volume, VolumeVector vector, ofstream &file) {
	auto nx = volume.nx;
	auto ny = volume.ny;
	auto nz = volume.nz;

	auto vx = vector.vx;
	auto vy = vector.vy;
	auto vz = vector.vz;

	file << format("VECTORS {} double", name) << endl;

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int p = j + i * ny;
				int id = k + p * nz;

				file << format("{:.{}f} {:.{}f} {:.{}f}", vx[id], precision, vy[id], precision, vz[id], precision) << endl;
			}
		}
	}
}

void Writer::writeVolumeScalar(string name, size_t precision, VolumeState volume, const vector<double> &data, ofstream &file) {
	auto nx = volume.nx;
	auto ny = volume.ny;
	auto nz = volume.nz;
	
	file << format("SCALARS {} double", name) << endl;
	file << "LOOKUP_TABLE default" << endl;

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int p = j + i * ny;
				int id = k + p * nz;

				file << format("{:.{}f}", data[id], precision) << endl;
			}
		}
	}
}

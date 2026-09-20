#include "utils.h"

using namespace Writers;

void Utils::defaultWriteSurfaceDataHeader(
	TimeState time,
	SurfaceState surface,
	ofstream file
) {
    auto nx = surface.nx;
    auto ny = surface.ny;

    auto dx = surface.dx;
    auto dy = surface.dy;

    file << "# vtk DataFile Version 3.0" << endl;
    file << format("TIME {:.3f}", time.t) << endl;
    file << "ASCII" << endl;
    file << "DATASET STRUCTURED_GRID" << endl;
    file << format("DIMENSIONS {} {} 1", nx, ny) << endl;
    file << format("POINTS {} double", nx * ny) << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = dx * i;
            double y = dy * j;

            file << format("{:.3f} {:.3f} 0.0", x, y) << endl;
        }
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 double" << endl;
    file << format("{:.3f}", time.t) << endl;
    file << format("POINT_DATA {}", nx * ny) << endl;
}

void Utils::writeVolumeTopHeader(TimeState time, VolumeState volume, ofstream& file) {
    auto nx = volume.nx;
    auto ny = volume.ny;
    auto nz = volume.nz;

    file << "# vtk DataFile Version 3.0" << endl;
    file << format("TIME {:.3f}", time.t) << endl;
    file << "ASCII" << endl;
    file << "DATASET STRUCTURED_GRID" << endl;
    file << format("DIMENSIONS {} {} {}", nx, ny, nz) << endl;
    file << format("POINTS {} double", nx * ny * nz) << endl;
}

void Utils::writeVolumeBottomHeader(TimeState time, VolumeState volume, ofstream& file) {
    auto nx = volume.nx;
    auto ny = volume.ny;
    auto nz = volume.nz;

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 double" << endl;
    file << format("{:.3f}", time.t) << endl;
    file << format("POINT_DATA {}", nx * ny * nz) << endl;
}

void Utils::writeVolumeElemHeader(TimeState time, VolumeState volume, ofstream& file) {
    auto nx = volume.nx;
    auto ny = volume.ny;
    auto nz = volume.nz;

    auto dx = volume.dx;
    auto dy = volume.dy;

    auto dz = volume.dz;
    auto h = volume.h;

    writeVolumeTopHeader(time, volume, file);

    vector<double> z(nx * ny, 0);

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double x = dx * i;
                double y = dy * j;

                int p = j + i * ny;

                z[p] += h[p] * dz[k] / 2;

                file << format("{:.3f} {:.3f} {:.3f}", x, y, z[p]) << endl;

                z[p] += h[p] * dz[k] / 2;
            }
        }
    }

    writeVolumeBottomHeader(time, volume, file);
}

void Utils::writeVolumeNodeHeader(TimeState time, VolumeState volume, ofstream& file) {
    auto nx = volume.nx;
    auto ny = volume.ny;
    auto nz = volume.nz;

    auto dx = volume.dx;
    auto dy = volume.dy;

    auto dz = volume.dz;
    auto h = volume.h;

    writeVolumeTopHeader(time, volume, file);

    vector<double> z(nx * ny, 0);

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double x = dx * i;
                double y = dy * j;

                int p = j + i * ny;

                file << format("{:.3f} {:.3f} {:.3f}", x, y, z[p]) << endl;

                if (k < nz - 1) {
                    z[p] += h[p] * dz[k];
                }
            }
        }
    }

    writeVolumeBottomHeader(time, volume, file);
}

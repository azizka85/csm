#ifndef WRITERS_UTILS_H
#define WRITERS_UTILS_H

#include <format>

#include <fstream>

#include <calc/data.h>

using namespace std;

using namespace Calc::Data;

namespace Writers::Utils {
	void defaultWriteSurfaceDataHeader(
		TimeState time,
		const SurfaceState& surface,
	ofstream &file
	);

	void writeVolumeTopHeader(TimeState time, const VolumeState& volume, ofstream& file);
	void writeVolumeBottomHeader(TimeState time, const VolumeState& volume, ofstream& file);

	void writeVolumeElemHeader(TimeState time, const VolumeState& volume, ofstream& file);
	void writeVolumeNodeHeader(TimeState time, const VolumeState& volume, ofstream& file);
}

#endif 
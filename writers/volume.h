#ifndef WRITERS_VOLUME_H
#define WRITERS_VOLUME_H

#include <span>
#include <vector>
#include <variant>

#include <fstream>

#include <filesystem>

#include <calc/data.h>

using namespace std;
using namespace std::filesystem;

using namespace Calc::Data;

namespace Writers::Volume {
	using WriteVolumeDataHeaderFn = void(*) (
		TimeState time,
		VolumeState volume,
		ofstream file
	);

	using DataKind = variant<
		VolumeVector,
		vector<double>
	>;

	class Writer {
		private:
			string fileName;
			WriteVolumeDataHeaderFn writeVolumeDataHeader;

		public:
			Writer(string fileName, WriteVolumeDataHeaderFn writeVolumeDataHeader);

			template<size_t N>
			void write(
				TimeState time,
				VolumeState volume,
				span<Column, N> columns,
				span<DataKind, N> data,
				path outDir
			);

			void writeVolumeVector(string name, size_t precision, VolumeState volume, VolumeVector vector, ofstream &file);
			void writeVolumeScalar(string name, size_t precision, VolumeState volume, vector<double> data, ofstream &file);
	};	
}

#endif 
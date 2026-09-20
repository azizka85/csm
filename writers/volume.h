#ifndef WRITERS_VOLUME_H
#define WRITERS_VOLUME_H

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
		ofstream &file
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
			
			path write(
				TimeState time,
				VolumeState volume,
				const vector<Column> &columns,
				const vector<DataKind> &data,
				path outDir
			);

			void writeVolumeVector(string name, size_t precision, VolumeState volume, VolumeVector vector, ofstream &file);
			void writeVolumeScalar(string name, size_t precision, VolumeState volume, const vector<double> &data, ofstream &file);
	};	
}

#endif 
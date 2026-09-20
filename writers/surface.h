#ifndef WRITERS_SURFACE_H
#define WRITERS_SURFACE_H

#include <vector>
#include <variant>

#include <fstream>

#include <filesystem>

#include <calc/data.h>

using namespace std;
using namespace std::filesystem;

using namespace Calc::Data;

namespace Writers::Surface {	
	using WriteSurfaceDataHeaderFn = void(*) (
		TimeState time,
		SurfaceState surface,
		ofstream &file
	);

	using DataKind = variant<
		SurfaceVector,
		vector<double>
	>;

	class Writer {
		private:
			string fileName;
			WriteSurfaceDataHeaderFn writeSurfaceDataHeader;

		public:
			Writer(string fileName, WriteSurfaceDataHeaderFn writeSurfaceDataHeader);

			path write(
				TimeState time,
				SurfaceState surface,
				const vector<Column>& columns,
				const vector<DataKind>& data,
				path outDir
			);

			void writeSurfaceVector(string name, size_t precision, SurfaceState surface, SurfaceVector vector, ofstream &file);
			void writeSurfaceScalar(string name, size_t precision, SurfaceState surface, const vector<double> &data, ofstream &file);
	};
}

#endif 
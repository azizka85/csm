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
		const SurfaceState& surface,
		ofstream& file
	);

	using DataKind = variant<
		SurfaceVector,
		Scalar
	>;

	class Writer {
		private:
			string fileName;
			WriteSurfaceDataHeaderFn writeSurfaceDataHeader;

		public:
			Writer(string fileName, WriteSurfaceDataHeaderFn writeSurfaceDataHeader);

			path write(
				TimeState time,
				const SurfaceState& surface,
				const vector<Column>& columns,
				const vector<DataKind>& data,
				path outDir
			);

			void writeSurfaceVector(string name, size_t precision, const SurfaceState &surface, const SurfaceVector &vector, ofstream &file);
			void writeSurfaceScalar(string name, size_t precision, const SurfaceState &surface, const Scalar &data, ofstream &file);
	};
}

#endif 
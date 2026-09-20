#ifndef WRITERS_STATISTICS_H
#define WRITERS_STATISTICS_H

#include <span>
#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace Writers::Statistics {
	struct Column {
		string name;
		size_t precision;
	};

	class Writer {
		public:			
			template <size_t N>
			path write(
				span<Column, N> columns, 
				vector<span<double, N>> data,
				path outDir
			);
	};
}

#endif 
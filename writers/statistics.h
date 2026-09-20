#ifndef WRITERS_STATISTICS_H
#define WRITERS_STATISTICS_H

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
			path write(
				const vector<Column> &columns,
				const vector<vector<double>> &data,
				path outDir
			);
	};
}

#endif 
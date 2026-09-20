#include <fstream>

#include <utils/checks/fs.h>

#include "statistics.h"

using namespace Utils::Checks;

using namespace Writers::Statistics;

template <size_t N>
path Writer::write(
	span<Column, N> columns,
	vector<span<double, N>> data,
	path outDir
) {
	FS::checkPath(outDir);

	auto filePath = outDir / path("convergence.csv");

	ofstream file(filePath);

	for (size_t i = 0; i < N; i++) {
		if (i > 0) {
			file << ", ";
		}

		file << columns[i].name;
	}

	file << endl;

	for (size_t i = 0; i < data.size(); i++) {
		for (size_t j = 0; j < N; j++) {
			if (j > 0) {
				file << ", ";
			}

			file << format("{:.{}f}", data[i][j], columns[j].precision);
		}

		file << endl;
	}
}
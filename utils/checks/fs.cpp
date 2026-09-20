#include <format>

#include <stdexcept>

#include "fs.h"

using namespace std;

using namespace Utils::Checks;

void FS::checkPathEmpty(path& targetPath) {
	if (targetPath.empty()) {
		throw runtime_error(
			format("path should not be empty")
		);
	}
}

void FS::checkPath(path& targetPath) {
	checkPathEmpty(targetPath);

	if (!exists(targetPath)) {
		throw runtime_error(
			format("path {} does not exist", targetPath.string())
		);
	}
}
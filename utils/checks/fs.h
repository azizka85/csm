#ifndef UTILS_CHECKS_FS_H
#define UTILS_CHECKS_FS_H

#include <filesystem>

using namespace std::filesystem;

namespace Utils::Checks::FS {
	void checkPathEmpty(path& targetPath);
	void checkPath(path& targetPath);
}

#endif
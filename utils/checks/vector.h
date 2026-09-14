#ifndef UTILS_CHECKS_VECTOR_H
#define UTILS_CHECKS_VECTOR_H

#include <span>

using namespace std;

namespace Utils::Checks::Vector {
	size_t check(span<const double> l, span<const double> c, span<const double> r);

	void checkU(size_t n, span<const double> u);
}

#endif
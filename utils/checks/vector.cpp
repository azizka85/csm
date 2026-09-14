#include <format>

#include <stdexcept>

#include "vector.h"

using namespace Utils::Checks;

size_t Vector::check(span<const double> l, span<const double> c, span<const double> r) {
	int n = c.size();

	if (l.size() != n) {
		throw runtime_error(
			format("The lengths of the l and c should be equal, but now the length of the l is {} and c is {}", l.size(), n)
		);
	}

	if (r.size() != n) {
		throw runtime_error(
			format("The lengths of the r and c should be equal, but now the length of the r is {} and c is {}", r.size(), n)
		);
	}
	
	return 0; 
}

void Vector::checkU(size_t n, span<const double> u) {
	if (u.size() != n) {
		throw runtime_error(
			format("The lengths of the u should be {}, but now the length of the u is {}", n, u.size())
		);
	}
}
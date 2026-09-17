#ifndef GENERATORS_UTILS_H
#define GENERATORS_UTILS_H 

#include <calc/data.h>

using namespace Calc;

namespace Generators::Utils {
	Data::CartesianBound toCartesian(Data::GeoBound bound);
	Data::GeoBound toGeo(Data::CartesianBound bound);
}

#endif 
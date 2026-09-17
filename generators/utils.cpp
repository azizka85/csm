#include <GeographicLib/UTMUPS.hpp>

#include "utils.h"

using namespace GeographicLib;

using namespace Generators::Utils;

Data::CartesianBound toCartesian(Data::GeoBound bound) {
	vector<Data::GeoPoint> corners = {
		Data::GeoPoint {
			.lat = bound.latMin,
			.lon = bound.lonMin
		},
		Data::GeoPoint {
			.lat = bound.latMin,
			.lon = bound.lonMax
		},
		Data::GeoPoint {
			.lat = bound.latMax,
			.lon = bound.lonMin
		},
		Data::GeoPoint {
			.lat = bound.latMax,
			.lon = bound.lonMax
		}
	};

	vector<double> east;
	vector<double> north;

	double x, y;
	
	int zone;
	bool northp;

	for (const auto& pt : corners) {
		UTMUPS::Forward(pt.lat, pt.lon, zone, northp, x, y, 39);

		east.push_back(x);
		north.push_back(y);
	}

	auto minmaxX = minmax_element(east.begin(), east.end());
	auto minmaxY = minmax_element(north.begin(), north.end());

	return Data::CartesianBound {
		.minX = *minmaxX.first,
		.minY = *minmaxY.first,

		.maxX = *minmaxX.second,
		.maxY = *minmaxY.second
	};
}

Data::GeoBound toGeo(Data::CartesianBound bound) {
	vector<Data::CartesianPoint> corners = {
		Data::CartesianPoint {
			.x = bound.minX,
			.y = bound.minY
		},
		Data::CartesianPoint {
			.x = bound.maxX,
			.y = bound.minY
		},
		Data::CartesianPoint {
			.x = bound.minX,
			.y = bound.maxY
		},
		Data::CartesianPoint {
			.x = bound.maxX,
			.y = bound.maxY
		}
	};

	vector<double> lats;
	vector<double> lons;

	double lat, lon;

	for (const auto& pt : corners) {
		UTMUPS::Reverse(39, true, pt.x, pt.y, lat, lon);

		lats.push_back(lat);
		lons.push_back(lon);
	}

	auto minmaxLat = minmax_element(lats.begin(), lats.end());
	auto minmaxLon = minmax_element(lons.begin(), lons.end());

	return Data::GeoBound {
		.latMin = *minmaxLat.first,
		.lonMin = *minmaxLon.first,

		.latMax = *minmaxLat.second,
		.lonMax = *minmaxLon.second
	};
}

#include <GeographicLib/UTMUPS.hpp>

#include <netcdf>

#include <boost/math/interpolators/bilinear_uniform.hpp>

#include <utils/checks/fs.h>
#include <utils/checks/area.h>
#include <utils/checks/bathymetry.h>

#include "../utils.h"

#include "gebco_generator.h"

using namespace GeographicLib;

using namespace netCDF;

using namespace boost::math::interpolators;

using namespace Utils;
using namespace Generators::Bathymetry;

GEBCOGenerator::GEBCOGenerator(GEBCOGeneratorParams params) {
	auto& bound = params.bound;

	Checks::Area::checkLatMinMax(bound.latMin, bound.latMax);
	Checks::Area::checkLonMinMax(bound.lonMin, bound.lonMax);

	Checks::Bathymetry::checkMinDepth(params.minDepth);

	Checks::FS::checkPath(params.filePath);

	this->params = move(params);
}

path GEBCOGenerator::addDirectory(path outDir) {
	auto& bound = params.bound;	

	return outDir / path(
		format(
			"bathymetry, GEBCO, lat={}-{}, lon={}-{}",
			bound.latMin, bound.latMax,
			bound.lonMin, bound.lonMax
		)
	);
}

vector<double> Generators::Bathymetry::GEBCOGenerator::generate(double nx, double ny) {
	Checks::Area::checkNumElemNodes(nx, ny);

    NcFile bathymetryFile(params.filePath.string(), NcFile::read);

    NcVar latVar = bathymetryFile.getVar("lat");
    size_t latSize = latVar.getDim(0).getSize();
    vector<double> lats(latSize);

    latVar.getVar(lats.data());

    double latStep = lats.size() > 1 ? lats[1] - lats[0] : 0;

    NcVar lonVar = bathymetryFile.getVar("lon");
    size_t lonSize = lonVar.getDim(0).getSize();
    vector<double> lons(lonSize);

    lonVar.getVar(lons.data());

    double lonStep = lons.size() > 1 ? lons[1] - lons[0] : 0;

    NcVar elevVar = bathymetryFile.getVar("elevation");
    vector<double> elevations(latSize * lonSize);

    elevVar.getVar(elevations.data());

    auto bound = Utils::toCartesian(params.bound);

    double dx = (bound.maxX - bound.minX) / (nx - 1);
    double dy = (bound.maxY - bound.minY) / (ny - 1);

    auto interpFunc = bilinear_uniform(
        move(elevations),
        latSize, lonSize,
        lonStep, latStep,
        lons.front(), lats.front()
    );

    vector<double> depths(nx * ny);

    double epsilon = 0.003;

    double lat, lon;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            double x = bound.minX + dx * i;
            double y = bound.minY + dy * j;

            UTMUPS::Reverse(39, true, x, y, lat, lon);            

            if (
                lat >= lats.front() + epsilon && lat <= lats.back() - epsilon &&
                lon >= lons.front() + epsilon && lon <= lons.back() - epsilon
            ) {
                double elevation = interpFunc(lon, lat);

                if (elevation < params.refDepth - params.minDepth) {
                    depths[p] = abs(elevation - params.refDepth);
                }
            }
        }
    }

    return depths;
}

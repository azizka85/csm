#include <span>

#include <stdexcept>

#include <GeographicLib/UTMUPS.hpp>

#include <netcdf>

#include <boost/math/interpolators/bilinear_uniform.hpp>

#include <utils/checks/area.h>
#include <utils/checks/phys.h>
#include <utils/checks/fs.h>

#include "../utils.h"

#include "ecmwf_generator.h"

using namespace GeographicLib;

using namespace netCDF;

using namespace boost::math::interpolators;

using namespace Generators::Wind;

ECMWFGenerator::ECMWFGenerator(ECMWFGeneratorParams params) {
	auto& bound = params.bound;

	::Utils::Checks::Area::checkLatMinMax(bound.latMin, bound.latMax);
	::Utils::Checks::Area::checkLonMinMax(bound.lonMin, bound.lonMax);	

	::Utils::Checks::FS::checkPath(params.filePath);

	this->params = move(params);
}

path Generators::Wind::ECMWFGenerator::addDirectory(path outDir) {
    auto attr = params.dirAttr.empty() ? "wind, ECMWF" : params.dirAttr;

	auto& bound = params.bound;

	return outDir / path(
        params.windStress->dirName(
            format(
                "{}, lat={}-{}, lon={}-{}", 
                attr, bound.latMin, bound.latMax, bound.lonMin, bound.lonMax
            )
        )
	);
}

vector<Data::Wind> Generators::Wind::ECMWFGenerator::generate(size_t nx, size_t ny) {
	::Utils::Checks::Area::checkNumElemNodes(nx, ny);

    NcFile windFile(params.filePath.string(), NcFile::read);

    NcVar timeVar = windFile.getVar("valid_time");
    size_t timeSize = timeVar.getDim(0).getSize();

    vector<int64_t> times(timeSize);

    timeVar.getVar(times.data());

    NcVar latVar = windFile.getVar("latitude");
    size_t latSize = latVar.getDim(0).getSize();
    vector<double> lats(latSize);

    latVar.getVar(lats.data());

    double latStep = lats.size() > 1 ? lats[0] - lats[1] : 0;

    NcVar lonVar = windFile.getVar("longitude");
    size_t lonSize = lonVar.getDim(0).getSize();
    vector<double> lons(lonSize);

    lonVar.getVar(lons.data());

    double lonStep = lons.size() > 1 ? lons[1] - lons[0] : 0;

    NcVar u10Var = windFile.getVar("u10");
    vector<double> u10Arr(timeSize * latSize * lonSize);

    u10Var.getVar(u10Arr.data());

    NcVar v10Var = windFile.getVar("v10");
    vector<double> v10Arr(timeSize * latSize * lonSize);

    v10Var.getVar(v10Arr.data());

    auto bound = Utils::toCartesian(params.bound);

    double dx = (bound.maxX - bound.minX) / (nx - 1);
    double dy = (bound.maxY - bound.minY) / (ny - 1);  

    vector<Data::Wind> res(timeSize);

    for (int k = 0; k < timeSize; k++) {
        double* u10SlicePtr = u10Arr.data() + k * latSize * lonSize;
        span<double> u10SliceSpan(u10SlicePtr, latSize * lonSize);

        auto u10Func = bilinear_uniform(
            move(u10SliceSpan),
            latSize, lonSize,
            lonStep, latStep,
            lons.front(), lats.back()
        );

        double* v10SlicePtr = v10Arr.data() + k * latSize * lonSize;
        span<double> v10SliceSpan(v10SlicePtr, latSize * lonSize);

        auto v10Func = bilinear_uniform(
            move(v10SliceSpan),
            latSize, lonSize,
            lonStep, latStep,
            lons.front(), lats.back()
        );

        vector<double> u10(nx * ny);
        vector<double> v10(nx * ny);

        vector<double> qx(nx * ny);
        vector<double> qy(nx * ny);

        int64_t refTime = timeSize > 0 ? times[0] : 0;

        double lat, lon;

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                int id = j + i * ny;

                double x = bound.minX + dx * i;
                double y = bound.minY + dy * j;

                UTMUPS::Reverse(39, true, x, y, lat, lon);

                double vu10 = 0;
                double vv10 = 0;

                double vqx = 0;
                double vqy = 0;

                if (
                    lat >= lats.back() && lat <= lats.front() &&
                    lon >= lons.front() && lon <= lons.back()

                ) {
                    vu10 = u10Func(lon, lat);
                    vv10 = v10Func(lon, lat);

                    auto stress = params.windStress->calculate(SpeedVector{ .u10 = vu10, .v10 = vv10 });

                    vqx = stress.qx;
                    vqy = stress.qy;
                }

                u10[id] = vu10;
                v10[id] = vv10;

                qx[id] = vqx;
                qy[id] = vqy;
            }
        }

        res[k].time = times[k] - refTime;

        res[k].speed.u10 = u10;
        res[k].speed.v10 = v10;

        res[k].stress.qx = qx;
        res[k].stress.qy = qy;
    }

    return res;
}

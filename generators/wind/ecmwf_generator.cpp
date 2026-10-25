#include <span>

#include <stdexcept>

#include <netcdf>

#include <boost/math/interpolators/bilinear_uniform.hpp>

#include <utils/checks/area.h>
#include <utils/checks/phys.h>
#include <utils/checks/fs.h>

#include "../utils.h"

#include "ecmwf_generator.h"

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

path ECMWFGenerator::addDirectory(path outDir) {
    auto attr = params.dirAttr.empty() ? "wind, ECMWF" : params.dirAttr;

	auto& bound = params.bound;

	return outDir / path(
        params.windStress->dirName(
            format(
                "{}", 
                attr
            )
        )
	);
}

vector<Data::Wind> ECMWFGenerator::generate(size_t nx, size_t ny, double dx, double dy) {
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

    auto flipRows = [&](const double* src) {
        vector<double> out(latSize * lonSize);

        for (size_t r = 0; r < latSize; r++) {
            copy(src + (latSize - 1 - r) * lonSize,
                src + (latSize - r) * lonSize,
                out.begin() + r * lonSize);
        }            

        return out;
    };

    auto bound = Utils::toCartesian(params.bound);

    vector<Data::Wind> res(timeSize);

    for (int k = 0; k < timeSize; k++) {
        auto u10Func = bilinear_uniform(
            flipRows(u10Arr.data() + k * latSize * lonSize),
            latSize, lonSize, lonStep, latStep, lons.front(), lats.back()
        );

        auto v10Func = bilinear_uniform(
            flipRows(v10Arr.data() + k * latSize * lonSize),
            latSize, lonSize, lonStep, latStep, lons.front(), lats.back()
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

                Utils::utm39Reverse(x, y, lat, lon);

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

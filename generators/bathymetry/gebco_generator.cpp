#include <queue>

#include <netcdf>

#include <boost/math/interpolators/bilinear_uniform.hpp>

#include <utils/checks/fs.h>
#include <utils/checks/area.h>
#include <utils/checks/bathymetry.h>

#include "../utils.h"

#include "gebco_generator.h"

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

	return outDir / path("bathymetry, GEBCO");
}

Data::Bathymetry GEBCOGenerator::generate(
    size_t nx, size_t ny,
    double dx, double dy
) {
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

    auto interpFunc = bilinear_uniform(
        move(elevations),
        latSize, lonSize,
        lonStep, latStep,
        lons.front(), lats.front()
    );

    Data::Bathymetry bathymetry{
        .h = vector<double>(nx * ny),
        .x = vector<double>(nx * ny),
        .y = vector<double>(nx * ny)
    };

    double epsilon = 0.003;

    const int ns = std::max(1, (int)std::ceil(std::max(dx, dy) / 1000.0));

    vector<double> meanDepth(nx * ny, 0.0);

    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            size_t p = j + i * ny;

            double x = bound.minX + dx * i;
            double y = bound.minY + dy * j;

            double lat, lon;
            Utils::utm39Reverse(x, y, lat, lon);

            bathymetry.x[p] = lon;
            bathymetry.y[p] = lat;

            double sum = 0;
            int cnt = 0;

            for (int a = 0; a < ns; a++) {
                for (int b = 0; b < ns; b++) {
                    double xs = x + dx * ((a + 0.5) / ns - 0.5);
                    double ys = y + dy * ((b + 0.5) / ns - 0.5);

                    double las, los;
                    Utils::utm39Reverse(xs, ys, las, los);

                    if (las >= lats.front() + epsilon && las <= lats.back() - epsilon &&
                        los >= lons.front() + epsilon && los <= lons.back() - epsilon) {
                        sum += std::max(params.refDepth - interpFunc(los, las), 0.0);
                        cnt++;
                    }
                }
            }

            if (cnt > 0) meanDepth[p] = sum / cnt;
        }
    }

    size_t seed = std::distance(meanDepth.begin(), std::max_element(meanDepth.begin(), meanDepth.end()));

    vector<char> connected(nx * ny, 0);
    std::queue<size_t> queue;

    if (meanDepth[seed] > params.minDepth) {
        connected[seed] = 1;
        queue.push(seed);
    }

    const long di[4] = { 1, -1, 0, 0 };
    const long dj[4] = { 0, 0, 1, -1 };

    while (!queue.empty()) {
        size_t p = queue.front();
        queue.pop();

        long i = (long)(p / ny), j = (long)(p % ny);

        for (int d = 0; d < 4; d++) {
            long ii = i + di[d], jj = j + dj[d];

            if (ii < 0 || jj < 0 || ii >= (long)nx || jj >= (long)ny) continue;

            size_t q = (size_t)jj + (size_t)ii * ny;

            if (!connected[q] && meanDepth[q] > params.minDepth) {
                connected[q] = 1;
                queue.push(q);
            }
        }
    }

    for (size_t p = 0; p < nx * ny; p++) {
        if (connected[p]) {
            bathymetry.h[p] = meanDepth[p];
        }
    }

    return bathymetry;
}

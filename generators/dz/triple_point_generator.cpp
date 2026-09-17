#include "triple_point_generator.h"

using namespace Generators::DZ;

TriplePointGenerator::TriplePointGenerator(TriplePointGeneratorParams params) {
	this->params = move(params);
}

path TriplePointGenerator::addDirectory(path outDir) {
	auto attr = params.dirAttr.empty() ? "tpdz" : params.dirAttr;

	return outDir / params.vertStep->dirName(attr);
}

vector<double> TriplePointGenerator::generate() {
	auto& vertStep = params.vertStep;
	
	vector<double> dz;

	double z = 0;
	double cdz = vertStep->calculate(z);

	dz.push_back(cdz);

    while (z < 1) {
        z += cdz;

        if (z >= 1) {
            if (z > 1) {
                z -= cdz;

                auto n = dz.size();

                auto dzd = 1 - z;

                if (n > 1) {
                    auto dzp = dz[n - 2];

                    dz[n - 2] = (dzd + dzp) / 2.;
                    dz[n - 1] = (dzd + dzp) / 2.;
                }
                else {
                    dz[n - 1] = dzd;
                }
            }

            break;
        }

        cdz = vertStep->calculate(z);
        dz.push_back(cdz);
    }

	return dz;
}

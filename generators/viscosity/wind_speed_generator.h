#ifndef GENERATORS_VISCOSITY_WIND_SPEED_GENERATOR_H
#define GENERATORS_VISCOSITY_WIND_SPEED_GENERATOR_H

#include <calc/turbulence.h>

#include "../viscosity.h"

namespace Generators::Viscosity {
	class WindSpeedGenerator : IGenerator {
	private:
		unique_ptr<Turbulence::IFrictionalInfluence> frictionalInfluence;
		unique_ptr<Turbulence::ITidalSurfaceViscosity> tidalSurfaceViscosity;
		unique_ptr<Turbulence::IViscosity> viscosity;

		string dirAttr;

	public:
		WindSpeedGenerator(
			unique_ptr<Turbulence::IFrictionalInfluence> frictionalInfluence,
			unique_ptr<Turbulence::ITidalSurfaceViscosity> tidalSurfaceViscosity,
			unique_ptr<Turbulence::IViscosity> viscosity,
			string dirAttr
		);

		path addDirectory(path outDir) override;
		vector<double> generate(Data::VolumeState volume, Data::WindCurrent windCurrent) override;
	};
}

#endif 
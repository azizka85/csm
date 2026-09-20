#ifndef GENERATORS_VISCOSITY_WIND_SPEED_GENERATOR_H
#define GENERATORS_VISCOSITY_WIND_SPEED_GENERATOR_H

#include <calc/turbulence.h>

#include "../viscosity.h"

namespace Generators::Viscosity {
	class WindSpeedGenerator : public IGenerator {
	private:
		shared_ptr<Turbulence::IFrictionalInfluence> frictionalInfluence;
		shared_ptr<Turbulence::ITidalSurfaceViscosity> tidalSurfaceViscosity;
		shared_ptr<Turbulence::IViscosity> viscosity;

		string dirAttr;

	public:
		WindSpeedGenerator(
			shared_ptr<Turbulence::IFrictionalInfluence> frictionalInfluence,
			shared_ptr<Turbulence::ITidalSurfaceViscosity> tidalSurfaceViscosity,
			shared_ptr<Turbulence::IViscosity> viscosity,
			string dirAttr
		);

		path addDirectory(path outDir) override;
		vector<double> generate(Data::VolumeState volume, Data::WindCurrent windCurrent) override;
	};
}

#endif 
#ifndef CONFIG_H
#define CONFIG_H

#include "configInterface.h"
#include "particleUtils.h"

using namespace std;
using namespace Particles;

namespace Config{

	const int DIM = 4;

	class Config : public ConfigInterface{
		public:
			Config() : ConfigInterface(){
				magneticField = "Tokamak";
				system = "GuidingCenter";
				integrator = "SymplecticExplicit1";

				h = 0.1;
				max_t = 1.E5;
				time_offset = 0;
				exit_on_error = false;
				error_threshold = 0.1;
				orbit_normalize = 50;

				z0.q << 0.050000000000000, 0.00000000000000 ,0.000000000000000 ,0.000390000000000;
				initialization_type = INIT_HAMILTONIAN;

				B0 = 1.;
				R0 = 1.;
				kt = 1.;
				q = 1.;

			};
			~Config(){};

			PositionMomentumPoint<DIM> z0;
	};
}

#endif
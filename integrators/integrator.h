#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "../systems/system.h"
#include "../systems/systemFactory.h"
#include <stdexcept>

using namespace Systems;

namespace Integrators{
	template <int DIM> class Integrator
	{
		public:
			Integrator(Config::Config* config);
			~Integrator(){};
			
			virtual PositionMomentumPoint<DIM> StepForward(PositionMomentumPoint<DIM> z) = 0;

			System<DIM>* system;

			double h;
	};

	template<int DIM> Integrator<DIM>::Integrator(Config::Config* config){
		if (config->h==0) throw invalid_argument("Invalid timestep h");
		system = systemFactory<DIM>(config->system,config);	

		h = config->h;		
	}
}

#endif
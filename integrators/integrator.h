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
			
			virtual PositionMomentumPoint<DIM> StepForward(PositionMomentumPoint<DIM> z, double h) = 0;

			System<DIM>* system;

			virtual PositionMomentumTwoPoints<DIM> initialize(PositionMomentumTwoPoints<DIM> z, initializationType init, double h);
	};

	template <int DIM> PositionMomentumTwoPoints<DIM> Integrator<DIM>::initialize(PositionMomentumTwoPoints<DIM> z, initializationType init, double h){
		if (init==INIT_MANUAL_POSITION_MOMENTUM) return z;
		else if (init==INIT_HAMILTONIAN){
			PositionPoints<DIM> q;
			q.q0 = z.q0;
			q.q1 = z.q1;

			z.q1 = z.q0;
			z.p1 = system->momentum(q);

			return z;
		}
		else throw invalid_argument("invalid initialization type");
	}

	template<int DIM> Integrator<DIM>::Integrator(Config::Config* config){
		if (config->h==0) throw invalid_argument("Invalid timestep h");
		system = systemFactory<DIM>(config->system,config);	
	}
}

#endif
#ifndef EXPLICIT_H
#define EXPLICIT_H

#include "integrator.h"
#include "../config.h"

namespace Integrators{
	template <int DIM> class ExplicitIntegrator: public Integrator<DIM>{
		public:
			ExplicitIntegrator(Config::Config* config): Integrator<DIM>::Integrator(config){};
			~ExplicitIntegrator(){};
	};
}

#endif
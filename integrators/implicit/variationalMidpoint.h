//build a variational integrator by discretizing the lagrangian with the midpoint rule
#ifndef VARIATIONALMIDPOINT_H
#define VARIATIONALMIDPOINT_H

#include "../variationalTrapezoidal.h"

namespace Integrators{
	template <int DIM> class VariationalMidpoint: public VariationalTrapezoidal<DIM>
	{
		public:
			VariationalMidpoint(Config::Config* config): VariationalTrapezoidal<DIM>(config, 0.5){};
			~VariationalMidpoint(){};
	};
}

#endif
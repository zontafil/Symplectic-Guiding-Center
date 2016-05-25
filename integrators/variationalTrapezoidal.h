// build a variational integrator starting from the continuous lagrangian and 
// discretizing it with the rectangle method (alpha is the position of the height of the rectangle )

#ifndef VARIATIONALTRAPEZOIDAL_H
#define VARIATIONALTRAPEZOIDAL_H

#include "variationalDiscreteLagrangian.h"
#include "../utils/particleUtils.h"

namespace Integrators{
	template <int DIM> class VariationalTrapezoidal: public VariationalDiscreteLagrangian<DIM>
	{
		private:
			double DiscreteLagrangian(PositionPoints<DIM> q, double h);
			double alpha;
		public:
			VariationalTrapezoidal(Config::Config* config, double a): VariationalDiscreteLagrangian<DIM>(config), alpha(a){};
			~VariationalTrapezoidal(){};
	};

	template <int DIM> double VariationalTrapezoidal<DIM>::DiscreteLagrangian(PositionPoints<DIM> q, double h){
		Matrix<double,DIM/2,1> qalpha = (1.-alpha)*q.q0 + alpha*q.q1;
		Matrix<double,DIM/2,1> dq = (q.q1-q.q0)/h;

		return (h * VariationalIntegrator<DIM>::system->Lagrangian(qalpha,dq));
	}
}

#endif
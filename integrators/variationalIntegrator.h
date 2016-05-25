//variational integrator:
// there must be an hamiltonian system.
// the integrator flow is defined by legendre left and right transforms (see paragraph 4.4.2)

// --> (q1,p1) = legendreRight(LegendreLeftInverse(q0,p0))

#ifndef VARIATIONAL_H
#define VARIATIONAL_H

#include "../utils/particleUtils.h"
#include "integrator.h"

using namespace Particles;

namespace Integrators{
	template <int DIM> class VariationalIntegrator : public Integrator<DIM>
	{
		public:
			VariationalIntegrator(Config::Config* config);
			~VariationalIntegrator(){};

			Matrix<double,DIM,1> StepForward(Matrix<double,DIM,1> z, double h);

			virtual PositionMomentumPoint<DIM> LegendreRight(PositionPoints<DIM> q, double h) = 0;
			virtual PositionPoints<DIM> LegendreLeftInverse(PositionMomentumPoint<DIM> z, double h) = 0;

			HamiltonianSystem<DIM>* system;

			virtual PhaseSpacePoints<DIM> initialize(PhaseSpacePoints<DIM> z, initializationType init, double h, Config::Config* config);
	
	};


	template <int DIM> VariationalIntegrator<DIM>::VariationalIntegrator(Config::Config* config): Integrator<DIM>(config){
		system = hamiltonianSystemFactory<DIM>(config->system,config);		
	}

	template <int DIM> PhaseSpacePoints<DIM> VariationalIntegrator<DIM>::initialize(PhaseSpacePoints<DIM> z, initializationType init, double h, Config::Config* config){
		if (init==INIT_MANUAL_MULTISTEP){

			//q0,q1 are provided. Compute p1 with LegendreRight

			PositionPoints<DIM> q;
			q.q0 = z.z0.head(DIM/2);
			q.q1 = z.z1.head(DIM/2);
			PositionMomentumPoint<DIM> z1 = LegendreRight(q, h);
			z.z1.tail(DIM/2) = z1.p;

			return z;
		}
		else if (init==INIT_HAMILTONIAN){

			//q0, q1 is provided. Compute p0 using continuous conjugate momentum

			PositionPoints<DIM> q;
			q.q0 = z.z0.head(DIM/2);
			q.q1 = z.z1.head(DIM/2);

			z.z1.head(DIM/2) = z.z0.head(DIM/2);
			z.z1.tail(DIM/2) = system->momentum(q);

			return z;
		}
		else return Integrator<DIM>::initialize(z,init, h, config);
	}

	template <int DIM> Matrix<double,DIM,1> VariationalIntegrator<DIM>::StepForward(Matrix<double,DIM,1> z0, double h){
		PositionMomentumPoint<DIM> z0pm, z1pm;
		Matrix<double,DIM,1> ret;
		z0pm.q = z0.head(DIM/2);
		z0pm.p = z0.tail(DIM/2);
		z1pm = LegendreRight(LegendreLeftInverse(z0pm, h), h);

		ret.head(DIM/2) = z1pm.q;
		ret.tail(DIM/2) = z1pm.p;

		return ret;
	}

}
#endif
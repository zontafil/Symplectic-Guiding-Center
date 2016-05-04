#ifndef RK4_H
#define RK4_H

#include "integrator.h"
#include "../config.h"
#include "../utils/particleUtils.h"

using namespace Particles;

namespace Integrators{
	template <int DIM> class RK4 : public Integrator<DIM>{
		public:
			RK4(Config::Config* config) : Integrator<DIM>::Integrator(config){
				system = Integrator<DIM>::system;
			};
			~RK4();

			PositionMomentumPoint<DIM> StepForward(PositionMomentumPoint<DIM> z, double h);

			System<DIM>* system;
	};

	template <int DIM> PositionMomentumPoint<DIM> RK4<DIM>::StepForward(PositionMomentumPoint<DIM> z, double h){
		Matrix<double,2*DIM,1> k1,k2,k3,k4,z0,z1;
		
		z0.head(DIM) = z.q;
		z0.tail(DIM) = z.p;
		
		k1 = system->f_eq_motion(z0);
		k2 = system->f_eq_motion(z0+0.5*h*k1);
		k3 = system->f_eq_motion(z0+0.5*h*k2);
		k4 = system->f_eq_motion(z0+h*k3);
		
		z1 = z0 + 1./6.* h *(k1+2.*k2+2.*k3+k4);
		
		PositionMomentumPoint<DIM> ret;
		ret.q = z1.head(DIM);
		ret.p = z1.tail(DIM);

		return ret;
	}
}

#endif
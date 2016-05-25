//RK4 for a generic system
// the system must have implemented eqs of motion.

#ifndef RK4_H
#define RK4_H

#include "../integrator.h"
#include "../../config.h"
#include "../../utils/particleUtils.h"

using namespace Particles;

namespace Integrators{
	template <int DIM> class RK4 : public Integrator<DIM>{
		public:
			RK4(Config::Config* config) : Integrator<DIM>(config){
				system = Integrator<DIM>::system;
			};
			~RK4(){};

			Matrix<double,DIM,1> StepForward(Matrix<double,DIM,1> z, double h);

			System<DIM>* system;
	};

	template <int DIM> Matrix<double,DIM,1> RK4<DIM>::StepForward(Matrix<double,DIM,1> z, double h){
		Matrix<double,DIM,1> k1,k2,k3,k4;

		k1 = system->f_eq_motion(z);
		k2 = system->f_eq_motion(z+0.5*h*k1);
		k3 = system->f_eq_motion(z+0.5*h*k2);
		k4 = system->f_eq_motion(z+h*k3);
		
		return (z + 1./6.* h *(k1+2.*k2+2.*k3+k4));
	}
}

#endif
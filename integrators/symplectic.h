#ifndef SYMPLECTIC_H
#define SYMPLECTIC_H

#include "../utils/particleUtils.h"

using namespace Particles;

namespace Integrators{
	template <int DIM> class SymplecticIntegrator : public Integrator<DIM>
	{
		public:
			SymplecticIntegrator(Config::Config* config): Integrator<DIM>(config){};
			~SymplecticIntegrator(){};

			PositionMomentumPoint<DIM> StepForward(PositionMomentumPoint<DIM> z, double h);

			virtual PositionMomentumPoint<DIM> LegendreRight(PositionPoints<DIM> q, double h) = 0;
			virtual PositionPoints<DIM> LegendreLeftInverse(PositionMomentumPoint<DIM> z, double h) = 0;

			PositionMomentumTwoPoints<DIM> initialize(PositionMomentumTwoPoints<DIM> z, initializationType init, double h);
	
	};

	template <int DIM> PositionMomentumTwoPoints<DIM> SymplecticIntegrator<DIM>::initialize(PositionMomentumTwoPoints<DIM> z, initializationType init, double h){
		if (init==INIT_MANUAL_MULTISTEP){
			PositionPoints<DIM> q;
			q.q0 = z.q0;
			q.q1 = z.q1;
			PositionMomentumPoint<DIM> z1 = LegendreRight(q, h);
			z.q1 = z1.q;
			z.p1 = z1.p;

			return z;
		}
		else return Integrator<DIM>::initialize(z,init, h);
	}

	template <int DIM> PositionMomentumPoint<DIM> SymplecticIntegrator<DIM>::StepForward(PositionMomentumPoint<DIM> z0, double h){
		return LegendreRight(LegendreLeftInverse(z0, h), h);
	}

}
#endif
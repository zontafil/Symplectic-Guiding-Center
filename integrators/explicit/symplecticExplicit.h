#ifndef SYMPLECTIC_H
#define SYMPLECTIC_H

#include "../../utils/particleUtils.h"
#include "../explicitIntegrator.h"

using namespace Particles;

namespace Integrators{
	template <int DIM> class SymplecticExplicitIntegrator : public ExplicitIntegrator<DIM>
	{
		public:
			SymplecticExplicitIntegrator(Config::Config* config): ExplicitIntegrator<DIM>(config){};
			~SymplecticExplicitIntegrator(){};

			PositionMomentumPoint<DIM> StepForward(PositionMomentumPoint<DIM> z, double h);

			virtual PositionMomentumPoint<DIM> LegendreRight(PositionPoints<DIM> q, double h) = 0;
			virtual PositionPoints<DIM> LegendreLeftInverse(PositionMomentumPoint<DIM> z, double h) = 0;

			PositionMomentumTwoPoints<DIM> initialize(PositionMomentumTwoPoints<DIM> z, initializationType init, double h);
	
	};

	template <int DIM> PositionMomentumTwoPoints<DIM> SymplecticExplicitIntegrator<DIM>::initialize(PositionMomentumTwoPoints<DIM> z, initializationType init, double h){
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

	template <int DIM> PositionMomentumPoint<DIM> SymplecticExplicitIntegrator<DIM>::StepForward(PositionMomentumPoint<DIM> z0, double h){
		return LegendreRight(LegendreLeftInverse(z0, h), h);
	}

}
#endif
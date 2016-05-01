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

			PositionMomentumPoint<DIM> StepForward(PositionMomentumPoint<DIM> z);

			virtual PositionMomentumPoint<DIM> LegendreRight(PositionPoints<DIM> q) = 0;
			virtual PositionPoints<DIM> LegendreLeftInverse(PositionMomentumPoint<DIM> z) = 0;
	
	};

	template <int DIM> PositionMomentumPoint<DIM> SymplecticIntegrator<DIM>::StepForward(PositionMomentumPoint<DIM> z0){
		return LegendreRight(LegendreLeftInverse(z0));
	}

}
#endif
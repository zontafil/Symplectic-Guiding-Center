#ifndef SYSTEM_H
#define SYSTEM_H
#include "../utils/particleUtils.h"

using namespace ParticleUtils;

namespace Systems{
	template <int DIM> class System{
		public:
			System();
			~System();

			PositionMomentumPoint<DIM> StepForward(PositionMomentumPoint<DIM> z);
			virtual PositionMomentumPoint<DIM> LegendreLeft(PositionPoints<DIM> q) = 0;
			virtual PositionMomentumPoint<DIM> LegendreRight(PositionPoints<DIM> q) = 0;
			virtual PositionPoints<DIM> LegendreLeftInverse(PositionMomentumPoint<DIM> z) = 0;
			virtual PositionPoints<DIM> LegendreRightInverse(PositionMomentumPoint<DIM> z) = 0;

			virtual double Hamiltonian(PositionMomentumPoint<DIM> z) = 0;

	};

	template <int DIM> System<DIM>::~System(){}
	template <int DIM> System<DIM>::System(){}

	template <int DIM> PositionMomentumPoint<DIM> System<DIM>::StepForward(PositionMomentumPoint<DIM> z0){
		return this->LegendreRight(this->LegendreLeftInverse(z0));
	}


}

#endif
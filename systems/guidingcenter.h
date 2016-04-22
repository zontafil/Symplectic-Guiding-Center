#ifndef GUIDINGCENTER_H
#define GUIDINGCENTER_H
#include "../utils/particleUtils.h"
#include "system.h"

using namespace ParticleUtils;

namespace Systems{

	template <int DIM> class GuidingCenter : public System<DIM>{ 
		public:
			GuidingCenter();
			~GuidingCenter();

			PositionMomentumPoint<DIM> LegendreLeft(PositionPoints<DIM> q);
			PositionMomentumPoint<DIM> LegendreRight(PositionPoints<DIM> q);
			PositionPoints<DIM> LegendreLeftInverse(PositionMomentumPoint<DIM> z);
			PositionPoints<DIM> LegendreRightInverse(PositionMomentumPoint<DIM> z);
	};

	template <int DIM> GuidingCenter<DIM>::GuidingCenter(){}
	template <int DIM> PositionMomentumPoint<DIM> GuidingCenter<DIM>::LegendreLeft(PositionPoints<DIM> q){}
	template <int DIM> PositionMomentumPoint<DIM> GuidingCenter<DIM>::LegendreRight(PositionPoints<DIM> q){}
	template <int DIM> PositionPoints<DIM> GuidingCenter<DIM>::LegendreLeftInverse(PositionMomentumPoint<DIM> z){}
	template <int DIM> PositionPoints<DIM> GuidingCenter<DIM>::LegendreRightInverse(PositionMomentumPoint<DIM> z){}

}

#endif
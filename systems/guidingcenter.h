#ifndef GUIDINGCENTER_H
#define GUIDINGCENTER_H
#include "../utils/particleUtils.h"
#include "system.h"

using namespace ParticleUtils;

namespace Systems{

	template <int DIM> struct EMfield{
		Matrix<double,DIM,1> B;
	};

	template <int DIM> class GuidingCenter : public System<DIM>{ 
		public:
			GuidingCenter();
			~GuidingCenter();

			PositionMomentumPoint<DIM> LegendreLeft(PositionPoints<DIM> q);
			PositionMomentumPoint<DIM> LegendreRight(PositionPoints<DIM> q);
			PositionPoints<DIM> LegendreLeftInverse(PositionMomentumPoint<DIM> z);
			PositionPoints<DIM> LegendreRightInverse(PositionMomentumPoint<DIM> z);

			double Hamiltonian(PositionMomentumPoint<DIM> z);

			// virtual EMfield<DIM> ComputeEMfield(Matrix<double,DIM,1> q) = 0;

		private:
			double mu;


	};

	template <int DIM> GuidingCenter<DIM>::GuidingCenter(){
		this->mu = 5;		
	}
	template <int DIM> PositionMomentumPoint<DIM> GuidingCenter<DIM>::LegendreLeft(PositionPoints<DIM> q){}
	template <int DIM> PositionMomentumPoint<DIM> GuidingCenter<DIM>::LegendreRight(PositionPoints<DIM> q){}
	template <int DIM> PositionPoints<DIM> GuidingCenter<DIM>::LegendreLeftInverse(PositionMomentumPoint<DIM> z){}
	template <int DIM> PositionPoints<DIM> GuidingCenter<DIM>::LegendreRightInverse(PositionMomentumPoint<DIM> z){}

	template <int DIM> double GuidingCenter<DIM>::Hamiltonian(PositionMomentumPoint<DIM> z){
		// Matrix<double,DIM,1> B = this->ComputeEMfield(z.q).B;
		// return (0.5*z.q(3)*z.q(3)+this->mu*B(z.q.head(3)).norm());
	}

}

#endif
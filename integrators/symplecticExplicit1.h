// Symplectic integrator for guiding center

#ifndef SYMPLECTICEXPLICIT1_H
#define SYMPLECTICEXPLICIT1_H

#include "symplectic.h"

using namespace ParticleUtils;

namespace Integrators{
	template <int DIM> class SymplecticExplicit1 : public SymplecticIntegrator<DIM>
	{
		public:
			SymplecticExplicit1(Config::Config* config);
			~SymplecticExplicit1(){};

			PositionMomentumPoint<DIM> LegendreLeft(PositionPoints<DIM> q);
			PositionMomentumPoint<DIM> LegendreRight(PositionPoints<DIM> q);
			PositionPoints<DIM> LegendreLeftInverse(PositionMomentumPoint<DIM> z);
			PositionPoints<DIM> LegendreRightInverse(PositionMomentumPoint<DIM> z);

			//override generic system with a guiding center system.
			//We need to use guiding center EM fields
			GuidingCenter<DIM>* system;

			double h;
			double mu;
	};

	template <int DIM> SymplecticExplicit1<DIM>::SymplecticExplicit1(Config::Config* config) : SymplecticIntegrator<DIM>(config){
		system = guidingcenterFactory<DIM>(config->system,config);		
		mu = system->mu;
		h = Integrator<DIM>::h;
	}

	template <int DIM> PositionMomentumPoint<DIM> SymplecticExplicit1<DIM>::LegendreLeft(PositionPoints<DIM> q){
			
	}
	template <int DIM> PositionMomentumPoint<DIM> SymplecticExplicit1<DIM>::LegendreRight(PositionPoints<DIM> q){

		PositionMomentumPoint<DIM> z;

		GuidingField field = system->fieldconfig->compute(q.q1);

		Matrix4d M;
		Vector4d W,Q;

		//~ //BUILD M
		M(0,1) = field.Bdag(2);
		M(0,2) = -field.Bdag(1);
		M(1,0) = -field.Bdag(2);
		M(1,2) = field.Bdag(0);
		M(2,0) = field.Bdag(1);  
		M(2,1) = -field.Bdag(0);   
		M(0,3) = -field.b(0);
		M(1,3) = -field.b(1);
		M(2,3) = -field.b(2);
		M(3,0)=field.b(0);
		M(3,1)=field.b(1);
		M(3,2)=field.b(2);
		M(0,0) = M(1,1) = M(2,2) = M(3,3) = 0;
		M/=2.;

		Vector4d dq = q.q1 - q.q0;
		Q = M*dq;

		z.q = q.q1;
		z.p.head(3) = Q.head(3) + field.Adag -h/2.*mu*field.B_grad;
		z.p(3) = Q(3) -h/2.*q.q1(3);

		return z;
	}
	template <int DIM> PositionPoints<DIM> SymplecticExplicit1<DIM>::LegendreLeftInverse(PositionMomentumPoint<DIM> z){

		PositionPoints<DIM> q;

		GuidingField field = system->fieldconfig->compute(z.q);
		  
		Matrix4d M;
		Vector4d W,Q;

		//~ //BUILD M
		M(0,1) = field.Bdag(2);
		M(0,2) = -field.Bdag(1);
		M(1,0) = -field.Bdag(2);
		M(1,2) = field.Bdag(0);
		M(2,0) = field.Bdag(1);  
		M(2,1) = -field.Bdag(0);   
		M(0,3) = -field.b(0);
		M(1,3) = -field.b(1);
		M(2,3) = -field.b(2);
		M(3,0)=field.b(0);
		M(3,1)=field.b(1);
		M(3,2)=field.b(2);
		M(0,0) = M(1,1) = M(2,2) = M(3,3) = 0;

		M/=2.;

		// //BUILD W
		W.head(3) = h/2.*(mu*field.B_grad) +field.Adag - z.p.head(3);
		W(3) = (h/2.*z.q(3) -z.p(3));

		Q = M.inverse()*W;

		q.q0 = z.q;
		q.q1.head(3) = (q.q0.head(3)+Q.head(3));
		q.q1(3) = (q.q0(3)+Q(3));

		return q;
	}
	template <int DIM> PositionPoints<DIM> SymplecticExplicit1<DIM>::LegendreRightInverse(PositionMomentumPoint<DIM> z){

	}

}
#endif
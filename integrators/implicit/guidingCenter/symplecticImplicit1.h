// Symplectic integrator for guiding center
// see paragraph 6.4.1

#ifndef SYMPLECTICIMPLICIT1_H
#define SYMPLECTICIMPLICIT1_H

#include "../../variationalImplicit.h"
#include <stdexcept>

using namespace Particles;

namespace Integrators{
	template <int DIM> class SymplecticImplicit1 : public VariationalImplicit<DIM>
	{
		public:
			SymplecticImplicit1(Config::Config* config);
			~SymplecticImplicit1(){};

			PositionMomentumPoint<DIM> LegendreRight(PositionPoints<DIM> q, double h);
			PositionMomentumPoint<DIM> LegendreLeft(PositionPoints<DIM> q, double h);

			//override generic system with a guiding center system.
			//We need to use guiding center EM fields
			GuidingCenter<DIM>* system;

			double mu;
	};

	template <int DIM> SymplecticImplicit1<DIM>::SymplecticImplicit1(Config::Config* config) : VariationalImplicit<DIM>(config){
		if (DIM!=8) throw invalid_argument("Invalid dimension for symplectic implicit 1: please use 8.");
		system = guidingcenterFactory<DIM>(config->system,config);		
		mu = config->mu;
	}

	template <int DIM> PositionMomentumPoint<DIM> SymplecticImplicit1<DIM>::LegendreRight(PositionPoints<DIM> q, double h){

		PositionMomentumPoint<DIM> z;

		Matrix<double,DIM/2,1> zm = ( q.q0 + q.q1) /2.;
		Vector3d xd = q.q1.head(3) - q.q0.head(3);
		double um = ( q.q0(3) + q.q1(3) ) /2.;

		GuidingField field = system->fieldconfig->compute(zm);

		z.q = q.q1;
		z.p.head(3) = 0.5*field.Adag_jac.transpose()*xd + field.Adag - h/2.*mu*field.B_grad;
		z.p(3) = 0.5*field.b.dot(xd) - h/2.*um;

		return z;
	}

	template <int DIM> PositionMomentumPoint<DIM> SymplecticImplicit1<DIM>::LegendreLeft(PositionPoints<DIM> q, double h){

		BOOST_LOG_TRIVIAL(debug) << "Implicit1_LegendreLeft: Legendre left:";
		BOOST_LOG_TRIVIAL(debug) << std::scientific << "Implicit1_LegendreLeft: q0: " << q.q0.transpose();
		BOOST_LOG_TRIVIAL(debug) << std::scientific << "Implicit1_LegendreLeft: q1: " << q.q1.transpose();

		PositionMomentumPoint<DIM> z;

		Matrix<double,DIM/2,1> qm = ( q.q0 + q.q1) /2.;
		Vector3d xd = q.q1.head(3) - q.q0.head(3);
		double um = ( q.q0(3) + q.q1(3) ) /2.;

		GuidingField field = system->fieldconfig->compute(qm);

		z.q = q.q0;
		z.p.head(3) = -0.5*field.Adag_jac.transpose()*xd + field.Adag + h/2.*mu*field.B_grad;
		z.p(3) = - 0.5*field.b.dot(xd) + h/2.*um;

		return z;
	}

}
#endif
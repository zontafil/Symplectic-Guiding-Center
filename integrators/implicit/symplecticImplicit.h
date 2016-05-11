#ifndef SYMPLECTICIMPLICIT_H
#define SYMPLECTICIMPLICIT_H

#include "../implicitIntegratorNewton.h"
#include "../../utils/particleUtils.h"

namespace Integrators{
	template <int DIM> class SymplecticImplicitIntegrator: public ImplicitIntegratorNewton<DIM>{
		private:
			//implicit function to invert
			Matrix<double,2*DIM,1> f_eq_motion_discrete(PositionMomentumTwoPoints<DIM> z, double h);

			virtual PositionMomentumPoint<DIM> LegendreRight(PositionPoints<DIM> q, double h) = 0;
			virtual PositionMomentumPoint<DIM> LegendreLeft(PositionPoints<DIM> q, double h) = 0;

		public:
			SymplecticImplicitIntegrator(Config::Config* config): ImplicitIntegratorNewton<DIM>::ImplicitIntegratorNewton(config){};
			~SymplecticImplicitIntegrator(){};	
	};

	template <int DIM> Matrix<double,2*DIM,1> SymplecticImplicitIntegrator<DIM>::f_eq_motion_discrete(PositionMomentumTwoPoints<DIM> z, double h){
		Matrix<double,2*DIM,1> f;
		PositionPoints<DIM> q;
		q.q0 = z.q0;
		q.q1 = z.q1;
		PositionMomentumPoint<DIM> z0 = LegendreLeft(q,h);
		PositionMomentumPoint<DIM> z1 = LegendreRight(q,h);

		f.head(DIM) = z.p0 - z0.p;
		f.tail(DIM) = z.p1 - z1.p;

		BOOST_LOG_TRIVIAL(debug) << "f_eq_motion_discrete(): f_p0_input:\t" << z.p0.transpose();
		BOOST_LOG_TRIVIAL(debug) << "f_eq_motion_discrete(): f_p_leg:\t" << z0.p.transpose();

		return f;
	}
}

#endif
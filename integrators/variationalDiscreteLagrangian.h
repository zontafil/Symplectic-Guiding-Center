//construct a variational integrator starting from a discrete lagrangian.
// the integrator is in general implicit, that's why it's an implementation of VariationalImplicit


#ifndef VARIATIONALIMPLICITDISCRETELAGRANGIAN_H
#define VARIATIONALIMPLICITDISCRETELAGRANGIAN_H

#include "variationalImplicit.h"

namespace Integrators{
	template <int DIM> class VariationalDiscreteLagrangian: public VariationalImplicit<DIM>
	{
		private:
			//discrete lagrangian
			virtual double DiscreteLagrangian(PositionPoints<DIM> q, double h) = 0;

			//space step for numerical derivatives
			const double hx;
		public:
			VariationalDiscreteLagrangian(Config::Config* config): VariationalImplicit<DIM>(config), hx(1.E-5){};
			~VariationalDiscreteLagrangian(){};

			//legendre left and right are computed numerically from discrete lagrangian.
			// see paragraph 4.2.2
			PositionMomentumPoint<DIM> LegendreLeft(PositionPoints<DIM> q, double h);
			PositionMomentumPoint<DIM> LegendreRight(PositionPoints<DIM> q, double h);
	};

	template <int DIM> PositionMomentumPoint<DIM> VariationalDiscreteLagrangian<DIM>::LegendreLeft(PositionPoints<DIM> q, double h){
		PositionMomentumPoint<DIM> ret;
		PositionPoints<DIM> dq0,dq1;
		ret.q = q.q0;

		for (unsigned int i = 0; i<DIM/2; i++){
			dq0 = dq1 = q;
			dq1.q0(i) += hx;
			dq0.q0(i) -= hx;

			ret.p(i) = 0.5 * ( DiscreteLagrangian(dq0,h) - DiscreteLagrangian(dq1,h) ) / hx;

		}		

		return ret;
	}

	template <int DIM> PositionMomentumPoint<DIM> VariationalDiscreteLagrangian<DIM>::LegendreRight(PositionPoints<DIM> q, double h){
		PositionMomentumPoint<DIM> ret;
		PositionPoints<DIM> dq0,dq1;
		ret.q = q.q1;

		for (unsigned int i = 0; i<DIM/2; i++){
			dq0 = dq1 = q;
			dq1.q1(i) += hx;
			dq0.q1(i) -= hx;

			ret.p(i) = 0.5 * ( DiscreteLagrangian(dq1,h) - DiscreteLagrangian(dq0,h) ) / hx;
		}		

		return ret;
	}

}

#endif
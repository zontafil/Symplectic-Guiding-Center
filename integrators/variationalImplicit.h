//variational integrator, implicit version.
//compute legendre left inverse using a first guess integrator and newton iterations

#ifndef VARIATIONALIMPLICIT_H
#define VARIATIONALIMPLICIT_H

#include "variationalIntegrator.h"
#include "../utils/particleUtils.h"
#include "explicitIntegratorFactory.h"

using namespace Particles;

namespace Integrators{
	template <int DIM> class VariationalImplicit: public VariationalIntegrator<DIM>
	{
		private:
			//first guess integrator for newton iterations
			Integrator<DIM>* firstGuess;

			//number of newton iterations
			unsigned int implicitIterations;

			PositionPoints<DIM> ImplicitIterationLegendreLeftInverse(PositionMomentumTwoPoints<DIM> z, double h);

			//space step for numerical derivatives
			const double hx;
		public:
			VariationalImplicit(Config::Config *config);
			~VariationalImplicit(){};

			virtual PositionMomentumPoint<DIM> LegendreLeft(PositionPoints<DIM> q, double h) = 0;

			PositionPoints<DIM> LegendreLeftInverse(PositionMomentumPoint<DIM> z, double h);
	};

	template <int DIM> VariationalImplicit<DIM>::VariationalImplicit(Config::Config* config): VariationalIntegrator<DIM>(config), hx(1.E-5){
		firstGuess = explicitIntegratorFactory<DIM>(config->first_guess_integrator,config);

		implicitIterations = config->implicit_iterations;
	}

	template <int DIM> PositionPoints<DIM> VariationalImplicit<DIM>::LegendreLeftInverse(PositionMomentumPoint<DIM> z0, double h){

		Matrix<double,DIM,1> z0vec,z1vec;
		z0vec.head(DIM/2) = z0.q;
		z0vec.tail(DIM/2) = z0.p;
		PositionMomentumTwoPoints<DIM> z;
		PositionMomentumPoint<DIM>z1;
		PositionPoints<DIM> q;

		//compute z1 using a first guess (i.e. RK4)
		z1vec = firstGuess->StepForward(z0vec,h);
		q.q0 = z0.q;
		q.q1 = z1vec.head(DIM/2);

		BOOST_LOG_TRIVIAL(debug) << std::scientific << "StepGuess:\t" << z1vec.transpose();

		for (unsigned int i = 0; i < implicitIterations; ++i)
		{
			z.q0 = z0.q;
			z.p0 = z0.p;
			z.q1 = q.q1;

			//compute q1 from q0,p0 and a first guess of q1
			q = ImplicitIterationLegendreLeftInverse(z,h);
		}

		return q;

	}

	template <int DIM> PositionPoints<DIM> VariationalImplicit<DIM>::ImplicitIterationLegendreLeftInverse(PositionMomentumTwoPoints<DIM> z, double h){
		//newton iteration for the legendre left inverse.
		// q1_new = q1 - f/f'

		Matrix<double,DIM/2,1> f, df1,df0, q1_new;
		Matrix<double,DIM/2,DIM/2> Jf;
		PositionPoints<DIM> q,dq1,dq0;
		q.q0 = z.q0;
		q.q1 = z.q1;

		PositionPoints<DIM> ret;
		ret.q0 = z.q0;

		f = z.p0 - LegendreLeft(q,h).p;

		BOOST_LOG_TRIVIAL(debug) << std::scientific << "------ ImplicitIteration() f:\t\t" << f.transpose();

		for (int j=0;j<DIM/2;j++){

			dq1 = dq0 = q;
			dq1.q1(j) += hx;
			dq0.q1(j) -= hx;

			df1 = z.p0 - LegendreLeft(dq1,h).p;
			df0 = z.p0 - LegendreLeft(dq0,h).p;

			Jf.col(j) = 0.5*(df1 - df0)/hx;
		}


		ret.q1 = z.q1 - Jf.inverse()*f;

		BOOST_LOG_TRIVIAL(debug) << "------ ImplicitIteration() q1_new:\t" << ret.q1.transpose();	

		return ret;
	}

}

#endif
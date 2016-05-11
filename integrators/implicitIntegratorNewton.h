#ifndef IMPLICITNEWTON_H
#define IMPLICITNEWTON_H

#include "integrator.h"
#include "../config.h"
#include "../utils/particleUtils.h"
#include "explicitIntegratorFactory.h"

namespace Integrators{
	template <int DIM> class ImplicitIntegratorNewton: public Integrator<DIM>{
		private:
			ExplicitIntegrator<DIM>* firstGuess;
			int implicitIterations;

			//implicit function to invert
			virtual Matrix<double,2*DIM,1> f_eq_motion_discrete(PositionMomentumTwoPoints<DIM> z, double h) = 0;

			PositionMomentumPoint<DIM> ImplicitIteration(PositionMomentumTwoPoints<DIM> z, double h);

			const double hx;
		public:
			ImplicitIntegratorNewton(Config::Config* config);
			~ImplicitIntegratorNewton(){};

			PositionMomentumPoint<DIM> StepForward(PositionMomentumPoint<DIM> z0, double h);
	};

	template <int DIM> ImplicitIntegratorNewton<DIM>::ImplicitIntegratorNewton(Config::Config* config): Integrator<DIM>::Integrator(config), hx(1.E-5){
		firstGuess = explicitIntegratorFactory<DIM>(config->first_guess_integrator,config);

		implicitIterations = config->implicit_iterations;
	}

	template <int DIM> PositionMomentumPoint<DIM> ImplicitIntegratorNewton<DIM>::StepForward(PositionMomentumPoint<DIM> z0, double h){
		PositionMomentumPoint<DIM> z1 = firstGuess->StepForward(z0,h);
		PositionMomentumTwoPoints<DIM> z;

		BOOST_LOG_TRIVIAL(debug) << std::scientific << "StepGuess:\t" << z1.q.transpose() << " " << z1.p.transpose();

		for (int i = 0; i < implicitIterations; ++i)
		{
			z.q0 = z0.q;
			z.p0 = z0.p;
			z.q1 = z1.q;
			z.p1 = z1.p;
			z1 = ImplicitIteration(z,h);
		}

		return z1;
	}

	template <int DIM> PositionMomentumPoint<DIM> ImplicitIntegratorNewton<DIM>::ImplicitIteration(PositionMomentumTwoPoints<DIM> z, double h){
		Matrix<double,2*DIM,1> z1, f, df1,df0, z_new;
		Matrix<double,2*DIM,2*DIM> Jf;
		PositionMomentumTwoPoints<DIM> dz1,dz0;
		PositionMomentumPoint<DIM> ret;

		z1.head(DIM) = z.q1;
		z1.tail(DIM) = z.p1;

		f = f_eq_motion_discrete(z,h);

		BOOST_LOG_TRIVIAL(debug) << std::scientific << "------ ImplicitIteration() f:\t\t" << f.transpose();

		for (int j=0;j<2*DIM;j++){

			dz1 = dz0 = z;
			if (j<DIM){
				dz1.q1(j) += hx;
				dz0.q1(j) -= hx;
			}
			else{
				dz1.p1(j-DIM) += hx;
				dz0.p1(j-DIM) -= hx;
			}

			df1 = f_eq_motion_discrete(dz1, h);
			df0 = f_eq_motion_discrete(dz0, h);

			Jf.col(j) = 0.5*(df1 - df0)/hx;
		}

		z_new = z1 - Jf.inverse()*f;
		ret.q = z_new.head(DIM);
		ret.p = z_new.tail(DIM);

		BOOST_LOG_TRIVIAL(debug) << "------ ImplicitIteration() z_new:\t" << z_new.transpose();	

		return ret;
	}
}

#endif
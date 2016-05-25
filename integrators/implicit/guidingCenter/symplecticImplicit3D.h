//Guiding center implicit integrator, scheme 2 (guiding center 3D lagrangian)
// see paragraph 6.4.2

#ifndef SYMPLECTICIMPLICIT2_H
#define SYMPLECTICIMPLICIT2_H

#include "../../variationalTrapezoidal.h"

namespace Integrators{
	template <int DIM> class SymplecticImplicit3D: public VariationalTrapezoidal<DIM>
	{
		private:
			Matrix<double,4,1> stepRK4(Matrix<double,4,1> z0, double h);
			Matrix<double,4,1> f_eq_motion(Matrix<double,4,1> z);
			const double mu;
		public:
			SymplecticImplicit3D(Config::Config* config): VariationalTrapezoidal<DIM>(config, 0.5), mu(config->mu){
				if (DIM!=6) throw invalid_argument("Invalid dimension for symplectic implicit 3D: please use 6.");
				system = guidingcenterFactory<DIM>(config->system,config);
			};
			~SymplecticImplicit3D(){};

			GuidingCenter<DIM>* system;

			PhaseSpacePoints<DIM> initialize(PhaseSpacePoints<DIM> z, initializationType init, double h, Config::Config* config);

	};

	template <int DIM> PhaseSpacePoints<DIM> SymplecticImplicit3D<DIM>::initialize(PhaseSpacePoints<DIM> z, initializationType init, double h, Config::Config* config){

		if (init==INIT_HAMILTONIAN){
			//compute p0 from x0 and u0 using continuous momenta equation (6.42)
			GuidingField field = system->fieldconfig->compute(z.z0.head(3));
			z.z1.head(3) = z.z0.head(3);
			z.z1.tail(3) = config->guidingcenter3D_u0*field.b + field.A;
			return z;
		}
		else if (init==INIT_LAGRANGIAN){

			// compute (x1,u1) from (x0,u0) with RK4 and p0 from (x0,x1) using legendre left

			Matrix<double,4,1> q1;
			q1.head(3) = z.z0.head(3);
			q1(3) = config->guidingcenter3D_u0;

			for (int i=0; i< config->init_steps; i++){
				q1 = stepRK4(q1,h/config->init_steps);
			}

			PositionPoints<DIM> q;
			q.q0 = z.z0.head(3);
			q.q1 = q1.head(3);

			PositionMomentumPoint<DIM> z0 = VariationalDiscreteLagrangian<DIM>::LegendreLeft(q,h);

			z.z1.head(3) = z0.q;
			z.z1.tail(3) = z0.p;

			return z;
		}
		else if (init==INIT_MANUAL_MULTISTEP){

			//compute (x0,p0) from (x0,x1) using Legendre left

			PositionPoints<DIM> q;
			q.q0 = z.z0.head(DIM/2);
			q.q1 = z.z1.head(DIM/2);
			PositionMomentumPoint<DIM> z1 = VariationalDiscreteLagrangian<DIM>::LegendreRight(q, h);
			z.z1.tail(DIM/2) = z1.p;

			return z;
		}
		else return Integrator<DIM>::initialize(z,init, h, config);



	}

	template <int DIM> Matrix<double,4,1> SymplecticImplicit3D<DIM>::stepRK4(Matrix<double,4,1> q0, double h){
		Matrix<double,4,1> k1,k2,k3,k4;

		k1 = f_eq_motion(q0);
		k2 = f_eq_motion(q0+0.5*h*k1);
		k3 = f_eq_motion(q0+0.5*h*k2);
		k4 = f_eq_motion(q0+h*k3);
		
		return (q0 + 1./6.* h *(k1+2.*k2+2.*k3+k4));

	}

	template <int DIM> Matrix<double,4,1> SymplecticImplicit3D<DIM>::f_eq_motion(Matrix<double,4,1> z){
		//equation of motion of guiding center 4D

		Matrix<double,4,1> f;

		GuidingField field = system->fieldconfig->compute(z.head(3));

		f.setZero();
		
		Matrix<double,3,1> E_dag;
		double B_dag_par;
		B_dag_par = field.Bdag.dot(field.b);
		E_dag = -mu*field.B_grad;
		
		f.head(3) = (z(3)* field.Bdag - field.b.cross(E_dag))/B_dag_par;
		f(3) = field.Bdag.dot(E_dag)/B_dag_par;

		return f;	
	}


}

#endif
//guiding center system, degenerate version
//		--> dimension must be 8 or 6
//		--> H(q,p) = H(q)
//		--> momentum(q0,q1) = momentum(q0)

// Use with variational phase-space symplectic integrators

#ifndef GUIDINGCENTER_H
#define GUIDINGCENTER_H
#include "../utils/particleUtils.h"
#include "../emfields/guidingfield.h"
#include "hamiltonianSystem.h"
#include <stdexcept>

using namespace Particles;
using namespace GuidingFields;

namespace Systems{

	template <int DIM> class GuidingCenter : public HamiltonianSystem<DIM>{ 
		private:
			static const double mu = 2.25E-6;
			static const double hx = 1.E-5;  //step for numerical derivative
		public:
			GuidingCenter(Config::Config* config);
			~GuidingCenter(){};

			double Hamiltonian(PositionMomentumPoint<DIM> z);
			double Lagrangian(Matrix<double,DIM/2,1> q, Matrix<double,DIM/2,1> v);
			Matrix<double,DIM/2,1> momentum(PositionPoints<DIM> q);
			Matrix<double,DIM,1> f_eq_motion(Matrix<double,DIM,1> z);

			GuidingFieldConfiguration<DIM> *fieldconfig;

	};


	template<int DIM> GuidingCenter<DIM>::GuidingCenter(Config::Config* config) : HamiltonianSystem<DIM>(config){
		if ((DIM!=8) && (DIM!=6)) throw invalid_argument("Invalid Guiding Center dimension: please use 8 or 6");

		//build an em guiding field.
		fieldconfig = new GuidingFieldConfiguration<DIM>(config);
	}

	template<int DIM> Matrix<double,DIM/2,1> GuidingCenter<DIM>::momentum(PositionPoints<DIM> q){

		//degenerate momenta for guiding center 8D.
		//see paragraph 6.3 of PDF

		Matrix<double,DIM/2,1> p;
		p.setZero();
		if (DIM==8){
			GuidingField field = fieldconfig->compute(q.q0);

			p(3) = 0.;
			p.head(3) = field.A+q.q0(3)*field.b;

		}
		return p;

	};
	template<int DIM> double GuidingCenter<DIM>::Hamiltonian(PositionMomentumPoint<DIM> z){
		if (DIM==8){
			//guiding center hamiltonian
			GuidingField field = this->fieldconfig->compute(z.q);
			return (0.5*z.q(3)*z.q(3)+mu*field.Bnorm);
		}
		else{
			//guiding center hamiltonian, 3D version. See parapraph 6.4.2
			GuidingField field = this->fieldconfig->compute(z.q);
			Matrix<double,3,1> p = z.p.head(3);
			double u = (p - field.A).norm();
			return (0.5*u*u+mu*field.Bnorm);
		}
	}
	template<int DIM> double GuidingCenter<DIM>::Lagrangian(Matrix<double,DIM/2,1> q, Matrix<double,DIM/2,1> v){
		if (DIM==8){
			//guiding center lagrangian
			GuidingField field = this->fieldconfig->compute(q);
			return (field.Adag.dot(v.head(3)) - (0.5*q(3)*q(3)+mu*field.Bnorm));
		}
		else{
			//guiding center lagrangian, 3D version. See parapraph 6.4.2
			Matrix<double,3,1> v3D = v.head(3);
			GuidingField field = this->fieldconfig->compute(q);
			double u = field.b.dot(v3D);
			return (field.A.dot(v3D) + 0.5*u*u - mu*field.Bnorm);
		}
	}

	template <int DIM> Matrix<double,DIM,1> GuidingCenter<DIM>::f_eq_motion(Matrix<double,DIM,1> z){
		Matrix<double,DIM,1> f;
		if (DIM==8){

			//guiding center eqs of motion. See paragraph 3.3

			GuidingField field = fieldconfig->compute(z.head(4));

			f.setZero();
			
			Matrix<double,3,1> E_dag;
			double B_dag_par;
			B_dag_par = field.Bdag.dot(field.b);
			E_dag = -mu*field.B_grad;
			
			f.head(3) = (z(3)* field.Bdag - field.b.cross(E_dag))/B_dag_par;
			f(3) = field.Bdag.dot(E_dag)/B_dag_par;

		}
		return f;	
	}
}

#endif
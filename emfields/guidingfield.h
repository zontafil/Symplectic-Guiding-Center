//Guiding Field Configuration.
// Expose useful guiding center quantities starting from a EM field
// i.e. gradient of B, A_dagger etc (see GuidingField struct)
#ifndef GUIDINGFIELD_H
#define GUIDINGFIELD_H

#include "emField.h"
#include "emFieldFactory.h"
#include <stdexcept>

using namespace EMFields;

namespace GuidingFields{

	//output type of computation
	struct GuidingField{
		Matrix<double,3,1> B, A, Adag, phi_grad, b, Bdag, B_grad;
		Matrix<double,3,3> Adag_jac, B_hessian, b_jac;
		double Bnorm, phi;
	};

	template <int DIM> class GuidingFieldConfiguration
	{
		private:
			Vector3d B_grad(Vector3d x);
			Matrix<double,3,3> B_hessian(Vector3d x);
			const double mu;
			const double hx;  //step for numerical derivative
			EMField* field;
		public:
			GuidingFieldConfiguration(Config::Config* config);
			~GuidingFieldConfiguration(){};

			//compute the field from q
			// q is (x,u) for 8D, (x) for 6D
			// if the dimension is 6, some quantities, like A_dagger are not computed
			GuidingField compute(Matrix<double,DIM/2,1> q);
	};

	template <int DIM> GuidingFieldConfiguration<DIM>::GuidingFieldConfiguration(Config::Config* config): mu(2.25E-6), hx(1.E-5){
		if ((DIM!=8) && (DIM!=6)) throw invalid_argument("Dimension must be 8 or 6 for the GuidingFieldConfiguration");

		field = EMFieldFactory(config->emField,config);
	};

	template <int DIM> Vector3d GuidingFieldConfiguration<DIM>::B_grad(Vector3d x){
		Vector3d dx(hx,hx,hx); //dx := (dx,dy,dz)
		Vector3d ret;
		Vector3d x0,x1;
		for (int j=0;j<3;j++){
		  x0 = x1 = x;
		  x0(j)-=dx(j);
		  x1(j)+=dx(j);
		  ret(j) = 0.5*(field->B(x1).norm() - field->B(x0).norm())/dx(j);
		}
		return ret;
	}
	
	template <int DIM> Matrix<double,3,3> GuidingFieldConfiguration<DIM>::B_hessian(Vector3d x){
		// COMPUTE B_HESSIAN
		Vector3d dx(hx,hx,hx); //dx := (dx,dy,dz)
		Matrix<double,3,3> ret;
		Vector3d x0,x1;
		for (int j=0;j<3;j++){
		  x0 = x1 = x;
		  x0(j)-=dx(j);
		  x1(j)+=dx(j);
		  ret.col(j) = 0.5*(B_grad(x1) - B_grad(x0))/dx(j);
		}

		return ret;
	}

	template <int DIM> GuidingField GuidingFieldConfiguration<DIM>::compute(Matrix<double,DIM/2,1> q){
		//TODO: implement phi

		BOOST_LOG_TRIVIAL(trace) << "Computing Magnetic field";
		BOOST_LOG_TRIVIAL(trace) << std::scientific << "q:\t\t" << q.transpose();
		
		GuidingField ret;
		Vector3d x = q.head(3);

		ret.A = field->A(x);
		ret.B = field->B(x);
		ret.Bnorm = ret.B.norm();
		ret.b = ret.B.normalized();

		double u;
		if (DIM==8) {
			u = q(3);
			ret.Adag = ret.A + u*ret.b;
		}
		else {
			ret.Adag.setZero();
			ret.Bdag.setZero();
			ret.Adag_jac.setZero();
		}

		Vector3d B0,x0,x1,B1;

		//COMPUTE GRADIENT(B),GRADIENT(phi),JAC(A_dagger)
		Matrix3d A_jac;
		for (int j=0;j<3;j++){
		  x0 = x1 = x;
		  x0(j)-=hx;
		  x1(j)+=hx;
		  B0 = field->B(x0);
		  B1 = field->B(x1);
		  
		  if (DIM==8) ret.Adag_jac.col(j) = 0.5*(field->A(x1) + u*B1.normalized() - field->A(x0)-u*B0.normalized())/hx;
		  A_jac.col(j) = 0.5*(field->A(x1) - field->A(x0))/hx;
		  ret.B_grad(j) = 0.5*(B1.norm() - B0.norm())/hx;
		  ret.b_jac.col(j) = 0.5*(B1.normalized() - B0.normalized())/hx;
		}

		//COMPUTE B_dagger
		if (DIM==8){
			ret.Bdag = ret.B;
			ret.Bdag(0) += u*(ret.b_jac(2,1)-ret.b_jac(1,2));
			ret.Bdag(1) += u*(ret.b_jac(0,2)-ret.b_jac(2,0));
			ret.Bdag(2) += u*(ret.b_jac(1,0)-ret.b_jac(0,1));		
		}

		ret.B_hessian = B_hessian(x);

		BOOST_LOG_TRIVIAL(trace) << std::scientific << "B:\t\t" << ret.B.transpose();
		BOOST_LOG_TRIVIAL(trace) << std::scientific << "A:\t\t" << ret.A.transpose();

		return ret;
	}
}

#endif
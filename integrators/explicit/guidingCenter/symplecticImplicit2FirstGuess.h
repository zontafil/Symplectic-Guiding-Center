// first guess for newton iteration of guiding center implicit 2 (3D version)

#ifndef IMPLICIT2FIRSTGUESS_H
#define IMPLICIT2FIRSTGUESS_H

namespace Integrators{
	template <int DIM> class SymplecticImplicit3DFirstGuess: public Integrator<DIM>
	{
		private:
			Matrix<double,4,1> f_eq_motion(Matrix<double,4,1> z);

			//override generic hamiltonian system with a guiding center system.
			//We need to use guiding center EM fields
			GuidingCenter<DIM>* system;
			const double mu;
		public:
			SymplecticImplicit3DFirstGuess(Config::Config* config): Integrator<DIM>(config), mu(config->mu){
				system = guidingcenterFactory<DIM>(config->system,config);		
			};
			~SymplecticImplicit3DFirstGuess(){};

			Matrix<double,DIM,1> StepForward(Matrix<double,DIM,1> z, double h);
	};

	template <int DIM> Matrix<double,DIM,1> SymplecticImplicit3DFirstGuess<DIM>::StepForward(Matrix<double,DIM,1> z, double h){

		// u0 = ||p0-A||/2 + mu*B
		// compute (x0,u0) --> (x1,u1) with RK4

		Matrix<double,4,1> k1,k2,k3,k4,q0,q1;

		q0.head(3) = z.head(3);
		GuidingField field = system->fieldconfig->compute(z.head(3));
		q0(3) = (z.tail(3) - field.A ).norm();

		k1 = f_eq_motion(q0);
		k2 = f_eq_motion(q0+0.5*h*k1);
		k3 = f_eq_motion(q0+0.5*h*k2);
		k4 = f_eq_motion(q0+h*k3);
		
		q1 = (q0 + 1./6.* h *(k1+2.*k2+2.*k3+k4));

		Matrix<double,DIM,1> ret;
		ret.head(3) = q1.head(3);
		return ret;
	}

	template <int DIM> Matrix<double,4,1> SymplecticImplicit3DFirstGuess<DIM>::f_eq_motion(Matrix<double,4,1> z){
		Matrix<double,4,1> f;

		GuidingField field = system->fieldconfig->compute(z.head(3));

		//compute Bdag manually, since u is not available in EMfields with dim=6
		double u = z(3);
		field.Bdag = field.B;
		field.Bdag(0) += u*(field.b_jac(2,1)-field.b_jac(1,2));
		field.Bdag(1) += u*(field.b_jac(0,2)-field.b_jac(2,0));
		field.Bdag(2) += u*(field.b_jac(1,0)-field.b_jac(0,1));		

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
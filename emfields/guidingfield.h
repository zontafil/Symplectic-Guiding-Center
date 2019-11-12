//Guiding Field Configuration.
// Expose useful guiding center quantities starting from a EM field
// i.e. gradient of B, A_dagger etc (see GuidingField struct)
#ifndef GUIDINGFIELD_H
#define GUIDINGFIELD_H

#include "emField.h"
#include "emFieldFactory.h"
#include <stdexcept>

namespace EMFields{

	//output type of computation
	struct GuidingField{
		Matrix<double,3,1> B, A, Adag, phi_grad, b, Bdag, B_grad;
		Matrix<double,3,3> Adag_jac, B_hessian, b_jac;
		double Bnorm, phi;
	};

	template <int DIM> class GuidingFieldConfiguration
	{
		public:
			GuidingFieldConfiguration(Config::Config* config);
			~GuidingFieldConfiguration(){};

			//compute the field from q
			virtual GuidingField compute(Matrix<double,DIM/2,1> q) = 0;
	};

	template <int DIM> GuidingFieldConfiguration<DIM>::GuidingFieldConfiguration(Config::Config* config) {
		if ((DIM!=8) && (DIM!=6)) throw invalid_argument("Dimension must be 8 or 6 for the GuidingFieldConfiguration");
	};

}

#endif
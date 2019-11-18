//Guiding Field Configuration.
// Expose useful guiding center quantities starting from a EM field
// i.e. gradient of B, A_dagger etc (see GuidingField struct)
#ifndef GUIDINGFIELD_H
#define GUIDINGFIELD_H

#include "ABfields/emField.h"
#include "ABfields/emFieldFactory.h"
#include <stdexcept>

namespace EMFields{

	//output type of computation
	struct GuidingField{
		Matrix<double,3,1> B, A, Adag, phi_grad, b, Bdag, B_grad, GradB_cyl, Bcyl;
		Matrix<double,3,3> Adag_jac, B_hessian, b_jac;
		double Bnorm, phi;
	};

	template <int DIM> class AB_dB_FieldBuilder
	{
		public:
			AB_dB_FieldBuilder(Config::Config* config);
			~AB_dB_FieldBuilder(){};

			//compute the field from q
			virtual GuidingField compute(Matrix<double,DIM/2,1> q) = 0;
	};

	template <int DIM> AB_dB_FieldBuilder<DIM>::AB_dB_FieldBuilder(Config::Config* config) {
		if ((DIM!=8) && (DIM!=6)) throw invalid_argument("Dimension must be 8 or 6 for the GuidingFieldConfiguration");
	};

}

#endif
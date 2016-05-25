// implicit guiding center version, modified by Qin. (paragraph 6.4.3)


#ifndef SEMIEXPLICITQIN_H
#define SEMIEXPLICITQIN_H

#include "../../variationalDiscreteLagrangian.h"

namespace Integrators{
	template <int DIM> class SemiexplicitQin: public VariationalDiscreteLagrangian<DIM>
	{
		private:
			double DiscreteLagrangian(PositionPoints<DIM> q, double h);
			const double mu;
		public:
			SemiexplicitQin(Config::Config* config): VariationalDiscreteLagrangian<DIM>(config), mu(config->mu){
				if (DIM!=8) throw invalid_argument("Invalid dimension for symplectic implicit 1: please use 8.");
				system = guidingcenterFactory<DIM>(config->system,config);		
			}
			~SemiexplicitQin(){};
			GuidingCenter<DIM>* system;

	};

	template <int DIM> double SemiexplicitQin<DIM>::DiscreteLagrangian(PositionPoints<DIM> q, double h){
		GuidingField field1 = system->fieldconfig->compute(q.q1);
		GuidingField field0 = system->fieldconfig->compute(q.q0);

		return (0.5*(field1.Adag+field0.Adag).dot(q.q1.head(3)-q.q0.head(3)) - 0.5*h*q.q1(3)*q.q0(3) - h*mu*field0.Bnorm  );
	}
}

#endif
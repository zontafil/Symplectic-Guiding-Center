//factory for integrators.
//use this if you want to compose an integrator inside a class (i.e. particle)

//i.e. if you want to create a RK4 integrator:
// Integrator<DIM> *rk = integratorFactory<DIM>("RK4",config);

#ifndef INTEGRATORFACTORY_H
#define INTEGRATORFACTORY_H

#include <stdexcept>
#include "integrator.h"
#include "explicit/guidingCenter/symplecticExplicit1.h"
#include "explicit/guidingCenter/symplecticExplicit2.h"
#include "explicit/guidingCenter/symplecticExplicit3.h"
#include "explicit/guidingCenter/symplecticExplicit4.h"
#include "explicit/RK4.h"

#include "implicit/guidingCenter/symplecticImplicit1.h"
#include "implicit/guidingCenter/symplecticSemiexplicitQin.h"
#include "implicit/guidingCenter/symplecticImplicit3D.h"
#include "implicit/variationalMidpoint.h"

using namespace std;

namespace Integrators{
	template <int DIM> Integrator<DIM> *integratorFactory(std::string const& integratorName, Config::Config* config){
		if (integratorName=="RK4") return new RK4<DIM>(config);
		else if (integratorName=="VariationalMidpoint") return new VariationalMidpoint<DIM>(config);
		else if (DIM==8){
			if (integratorName=="SymplecticExplicit1") return new SymplecticExplicit1<DIM>(config);
			else if (integratorName=="SymplecticExplicit2") return new SymplecticExplicit2<DIM>(config);
			else if (integratorName=="SymplecticExplicit3") return new SymplecticExplicit3<DIM>(config);
			else if (integratorName=="SymplecticExplicit4") return new SymplecticExplicit4<DIM>(config);
			else if (integratorName=="SymplecticImplicit1") return new SymplecticImplicit1<DIM>(config);
			else if (integratorName=="symplecticSemiexplicitQin") return new SemiexplicitQin<DIM>(config);
		}
		else if (DIM==6){
			if (integratorName=="SymplecticImplicit3D") return new SymplecticImplicit3D<DIM>(config);
		}

		throw invalid_argument("Invalid integrator "+ integratorName);
	}

}

#endif
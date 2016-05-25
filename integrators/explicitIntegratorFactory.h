//factory for explicit integrators.
//use this if you want to compose an explicit integrator inside a class (i.e. particle)

//i.e. if you want to create a RK4 integrator:
// Integrator<DIM> *rk = integratorFactory<DIM>("RK4",config);

//you want to use this instead of integrator factory, when you can't use an implicit integrator.
//i.e. a first guess newton integrator should be explicit.

#ifndef EXPLICITINTEGRATORFACTORY_H
#define EXPLICITINTEGRATORFACTORY_H

#include <stdexcept>
#include "integrator.h"
#include "explicit/guidingCenter/symplecticExplicit1.h"
#include "explicit/guidingCenter/symplecticExplicit2.h"
#include "explicit/guidingCenter/symplecticExplicit3.h"
#include "explicit/guidingCenter/symplecticExplicit4.h"
#include "explicit/guidingCenter/symplecticImplicit2FirstGuess.h"
#include "explicit/RK4.h"

using namespace std;

namespace Integrators{

	template <int DIM> Integrator<DIM> *explicitIntegratorFactory(std::string const& integratorName, Config::Config* config){
		if (integratorName=="RK4") return new RK4<DIM>(config);
		else if (DIM==8){
			if (integratorName=="SymplecticExplicit1") return new SymplecticExplicit1<DIM>(config);
			else if (integratorName=="SymplecticExplicit2") return new SymplecticExplicit2<DIM>(config);
			else if (integratorName=="SymplecticExplicit3") return new SymplecticExplicit3<DIM>(config);
			else if (integratorName=="SymplecticExplicit4") return new SymplecticExplicit4<DIM>(config);
		}
		else if (DIM==6){
			if (integratorName=="SymplecticImplicit3DFirstGuess") return new SymplecticImplicit3DFirstGuess<DIM>(config);
		}
	
		throw invalid_argument("Invalid explicit integrator "+ integratorName);
	}

}

#endif
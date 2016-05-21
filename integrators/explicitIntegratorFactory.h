#ifndef EXPLICITINTEGRATORFACTORY_H
#define EXPLICITINTEGRATORFACTORY_H

#include <stdexcept>
#include "integrator.h"
#include "explicit/guidingCenter/symplecticExplicit1.h"
#include "explicit/guidingCenter/symplecticExplicit2.h"
#include "explicit/guidingCenter/symplecticExplicit3.h"
#include "explicit/guidingCenter/symplecticExplicit4.h"
#include "explicit/RK4.h"

using namespace std;

namespace Integrators{

	template <int DIM> Integrator<DIM> *explicitIntegratorFactory(std::string const& integratorName, Config::Config* config){
		if (integratorName=="SymplecticExplicit1") return new SymplecticExplicit1<DIM>(config);
		else if (integratorName=="SymplecticExplicit2") return new SymplecticExplicit2<DIM>(config);
		else if (integratorName=="SymplecticExplicit3") return new SymplecticExplicit3<DIM>(config);
		else if (integratorName=="SymplecticExplicit4") return new SymplecticExplicit4<DIM>(config);
		else if (integratorName=="RK4") return new RK4<DIM>(config);
	
		throw invalid_argument("Invalid explicit integrator "+ integratorName);
	}

}

#endif
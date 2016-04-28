#ifndef INTEGRATORFACTORY_H
#define INTEGRATORFACTORY_H

#include <stdexcept>
#include "integrator.h"
#include "symplecticExplicit1.h"

using namespace std;

namespace Integrators{
	template <int DIM> Integrator<DIM> *integratorFactory(std::string const& integratorName, Config::Config* config){
		if (integratorName=="SymplecticExplicit1") return new SymplecticExplicit1<DIM>(config);
		else throw invalid_argument("Invalid integrator "+ integratorName);
	}

}

#endif
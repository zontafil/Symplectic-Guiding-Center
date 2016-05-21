#ifndef SYSTEMFACTORY_H
#define SYSTEMFACTORY_H

#include <stdexcept>
#include "system.h"
#include "guidingcenter.h"

using namespace std;

namespace Systems{
	
	template <int DIM> System<DIM> *systemFactory(std::string const& systemName, Config::Config* config){
		if (systemName=="GuidingCenter") return new GuidingCenter<DIM>(config);
		else throw invalid_argument("Invalid system "+ systemName);
	}
	template <int DIM> GuidingCenter<DIM> *guidingcenterFactory(std::string const& systemName, Config::Config* config){
		if (systemName=="GuidingCenter") return new GuidingCenter<DIM>(config);
		else throw invalid_argument("Invalid Guiding center system "+ systemName);
	}
	template <int DIM> HamiltonianSystem<DIM> *hamiltonianSystemFactory(std::string const& systemName, Config::Config* config){
		if (systemName=="GuidingCenter") return new GuidingCenter<DIM>(config);
		else throw invalid_argument("Invalid Hamiltonian system "+ systemName);
	}

}

#endif
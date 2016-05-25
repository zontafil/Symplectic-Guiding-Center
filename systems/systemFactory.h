//factory for systems.
//use this if you want to compose a system inside a class (i.e. particle)

//i.e. if you want to create a GuidingCenter system:
// System<DIM> *gc = systemFactory<DIM>("GuidingCenter",config);

#ifndef SYSTEMFACTORY_H
#define SYSTEMFACTORY_H

#include <stdexcept>
#include "system.h"
#include "guidingcenter.h"
#include "hamiltonianSystem.h"

using namespace std;

namespace Systems{
	
	template <int DIM> System<DIM> *systemFactory(std::string const& systemName, Config::Config* config){
		if (systemName=="GuidingCenter") return new GuidingCenter<DIM>(config);
		throw invalid_argument("Invalid system "+ systemName);
	}
	template <int DIM> GuidingCenter<DIM> *guidingcenterFactory(std::string const& systemName, Config::Config* config){
		if (systemName=="GuidingCenter") return new GuidingCenter<DIM>(config);
		throw invalid_argument("Invalid Guiding center system "+ systemName);
	}
	template <int DIM> HamiltonianSystem<DIM> *hamiltonianSystemFactory(std::string const& systemName, Config::Config* config){
		if (systemName=="GuidingCenter") return new GuidingCenter<DIM>(config);
		throw invalid_argument("Invalid Hamiltonian system "+ systemName);
	}

}

#endif
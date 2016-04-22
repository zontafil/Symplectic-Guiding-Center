#ifndef SYSTEMFACTORY_H
#define SYSTEMFACTORY_H

#include "system.h"
#include "guidingcenter.h"

namespace Systems{
	template <int DIM> System<DIM> *systemFactory(std::string const& systemName){
		System<DIM> *system = new GuidingCenter<DIM>();
		// if (systemName=="GuidingCenter") system = GuidingCenter<DIM>();

		return system;
	}
	
}

#endif
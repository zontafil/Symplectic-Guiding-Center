#ifndef SYSTEMFACTORY_H
#define SYSTEMFACTORY_H

#include <stdexcept>
#include "system.h"
#include "guidingcenter.h"

using namespace std;

namespace Systems{
	template <int DIM> System<DIM> *systemFactory(std::string const& systemName){
		if (systemName=="GuidingCenter") return new GuidingCenter<DIM>();
		else throw invalid_argument("Invalid system "+ systemName);
	}

}

#endif
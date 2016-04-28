#ifndef GUIDINGFIELD_FACTORY_H
#define GUIDINGFIELD_FACTORY_H

#include <stdexcept>
#include "tokamak.h"
#include "guidingfield.h"

using namespace std;

namespace GuidingFields{
	GuidingFieldConfiguration *guidingfieldFactory(std::string const& guidingfieldName, Config::Config* config){
		if (guidingfieldName=="Tokamak") return new Tokamak(config);
		else throw invalid_argument("Invalid guiding field "+ guidingfieldName);
	}

}

#endif
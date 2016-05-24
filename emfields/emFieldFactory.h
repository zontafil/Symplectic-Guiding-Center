#ifndef EMFIELDFACTORY_H
#define EMFIELDFACTORY_H

#include <stdexcept>
#include "emField.h"
#include "tokamak.h"
#include "forcefree.h"
#include "twoDimField.h"

using namespace std;

namespace EMFields{
	EMField *EMFieldFactory(std::string const& fieldName, Config::Config* config){
		if (fieldName=="Tokamak") return new Tokamak(config);
		else if (fieldName=="ForceFree") return new ForceFree(config);
		else if (fieldName=="TwoDimField") return new TwoDimField(config);

		throw invalid_argument("Invalid EMField "+ fieldName);
	}

}

#endif
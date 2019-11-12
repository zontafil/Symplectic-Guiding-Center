#ifndef GUIDINGFIELD_FACTORY_H
#define GUIDINGFIELD_FACTORY_H

#include <stdexcept>
#include "guidingfield.h"
#include "finiteDFromAB.h"
#include "splineField.h"

using namespace std;

namespace EMFields{
	template <int DIM> GuidingFieldConfiguration<DIM> *GuidingFieldFactory(std::string const& fieldName, Config::Config* config){
		if (fieldName=="finiteDFromAB") return new FiniteDFromAB<DIM>(config);
		else if (fieldName=="splineField") return new SplineField<DIM>(config);

		throw invalid_argument("Invalid Guiding Field Algorithm (config->guidingFieldAlgorithm). Choices: finiteDFromAB");
	}
}

#endif
#ifndef GUIDINGFIELD_FACTORY_H
#define GUIDINGFIELD_FACTORY_H

#include <stdexcept>
#include "AB_dB_Field.h"
#include "AB_dBfields/finiteDFromAB.h"
#include "AB_dBfields/finiteDFromA.h"
#include "AB_dBfields/splineField_BdB.h"

using namespace std;

namespace EMFields{
	template <int DIM> AB_dB_FieldBuilder<DIM> *AB_dB_FieldFactory(std::string const& fieldName, Config::Config* config){
		if (fieldName=="finiteDFromAB") return new FiniteDFromAB<DIM>(config);
		if (fieldName=="finiteDFromA") return new FiniteDFromA<DIM>(config);
		else if (fieldName=="splineField") return new SplineField_BdB<DIM>(config);

		throw invalid_argument("Invalid Guiding Field Algorithm (config->AB_dB_Algorithm). Choices: finiteDFromAB");
	}
}

#endif
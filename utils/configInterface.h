#ifndef CONFIGINTERFACE_H
#define CONFIGINTERFACE_H

#include <stdlib.h>
#include "particleUtils.h"

using namespace std;
using namespace Particles;

namespace Config{

	class ConfigInterface
	{
		public:
			string magneticField, system, integrator;
			double h;

			int max_t, time_offset, orbit_normalize;
			bool exit_on_error;
			double error_threshold;

			int print_precision, print_timestep_mult, file_timestep_mult, print_timestep_offset;

			string outFile;

			initializationType initialization_type;

			ConfigInterface(): h(0) {};
			~ConfigInterface(){};
	};

}

#endif
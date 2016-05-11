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
			string magneticField, system, integrator, auxiliary_integrator;
			double h;
			int init_steps;

			int max_t, time_offset, orbit_normalize;
			bool exit_on_error;
			double error_threshold;

			int print_precision, print_timestep_mult, file_timestep_mult, print_timestep_offset;

			string outFile;

			initializationType initialization_type;

			//force free specific
			bool forcefree_pert;
			double forcefree_kff; //perturbation magnitude
			double forcefree_a;

			string debug_level;

			ConfigInterface(): h(0), debug_level("info") {};
			~ConfigInterface(){};
	};

}

#endif
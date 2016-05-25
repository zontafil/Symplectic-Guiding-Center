//interface for configuration class
//see config.h for the implementation or copy a configuration file from "examples.h"

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
			//system to integrate (i.e. GuidingCenter)
			//see "systems/systemFactory.h" for a full list
			string system;

			//integrator to use (i.e. RK4, SymplecticExplicit3)
			//see "integrators/integratorFactory.h" for a full list
			string integrator;

			//timestep
			double h;

			// timestep to compute, initial t, number of steps per orbit
			int max_t, time_offset, orbit_normalize;

			//initialization type ({INIT_HAMILTONIAN, INIT_LAGRANGIAN, INIT_MANUAL, INIT_MANUAL_MULTISTEP})
			initializationType initialization_type;

			//auxiliary integrator and init steps for INIT_LAGRANGIAN (integrator to compute z1 from z0)
			string auxiliary_integrator;
			int init_steps;

			// first guess integrator for implicit and number of newton iterations (must be explicit)
			string first_guess_integrator;
			int implicit_iterations;
			
			//exit if the error is > error_threshold
			bool exit_on_error;
			double error_threshold;

			// print digit precision, print to screen and file every n steps
			int print_precision, print_timestep_mult, file_timestep_mult, print_timestep_offset;

			//output file
			string outFile;

			//info(default)/debug/trace
			string debug_level;

			//GUIDING CENTER SPECIFIC

				//EM field.
				//See "emfields/emFieldFactory.h" for a full list.
				string emField;

				double mu;

				//u0 for guiding center implicit 2 (3D)
				double guidingcenter3D_u0;

				//force free specific
				bool forcefree_pert;
				double forcefree_kff; //perturbation magnitude
				double forcefree_a;

			ConfigInterface(): h(0), debug_level("info") {};
			~ConfigInterface(){};
	};

}

#endif
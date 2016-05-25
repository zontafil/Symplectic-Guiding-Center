//configuraion file. See "utils/configInterface.h" for the documentation

#ifndef CONFIG_H
#define CONFIG_H

#include "utils/configInterface.h"
#include "utils/particleUtils.h"

using namespace std;
using namespace Particles;

namespace Config{

	//dimension of the phase space of the system
	const int DIM = 8;

	class Config : public ConfigInterface{
		public:
			Config() : ConfigInterface(){
				emField = "Tokamak"; //Tokamak,ForceFree,TwoDimField

				system = "GuidingCenter";
				integrator = "RK4";

				h = 800; //timestep
				max_t = 1E4; //number of timesteps to compute
				time_offset = 0; //timestep when to begin simulation
				exit_on_error = true; //exit if the error is too high
				error_threshold = 0.1; //error threshold to exit simulation if above
				orbit_normalize = 50; //timestep normalization

				z.z0.head(DIM/2) << 0.050000000000000, 0.00000000000000 ,0.000000000000000 ,0.000390000000000; //initial condition
				initialization_type = INIT_MANUAL; // {INIT_HAMILTONIAN, INIT_LAGRANGIAN, INIT_MANUAL, INIT_MANUAL_MULTISTEP}
				auxiliary_integrator = "RK4";
				init_steps = 1000;

				B0 = 1.; // magnetic field B0
				R0 = 1.; // major radius
				kt = 1.; // toroidal field correction
				q = 2.; // safety factor

				mu = 2.25E-6;

				outFile = "out/out.txt"; //output file

				print_precision = 10; // output digit precision
				print_timestep_mult = 10000; // print to screen every n timesteps
				file_timestep_mult = 1; //print to file every n timesteps
				print_timestep_offset = 1; //print to screen initial offset

				debug_level = "info"; //set to info(default)/debug/trace

			};
			~Config(){};

			PhaseSpacePoints<DIM> z;

			//guiding center specific
			double B0,R0,kt,q;
	};
}

#endif
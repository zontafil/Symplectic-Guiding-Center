#ifndef CONFIG_H
#define CONFIG_H

#include <stdlib.h>

using namespace std;

namespace Config{

	const int DIM = 4;

	class ConfigInterface
	{
		public:
			int DIM;
			string magneticField, system, integrator;
			double h;

			double B0,R0,kt,q;

			int max_t, time_offset, orbit_normalize;
			bool exit_on_error;
			double error_threshold;

			ConfigInterface(): DIM(Config::DIM), h(0) {};
			~ConfigInterface(){};
	};


	class Config : public ConfigInterface
	{
	public:
		Config() : ConfigInterface(){
			magneticField = "Tokamak";
			system = "GuidingCenter";
			integrator = "SymplecticExplicit1";
			h = 0.1;

			B0 = 1.;
			R0 = 1.;
			kt = 1.;
			q = 1.;

			max_t = 1.E5;
			time_offset = 0;
			exit_on_error = false;
			error_threshold = 0.1;
			orbit_normalize = 50;

			q0 << 0.050000000000000, 0.00000000000000 ,0.000000000000000 ,0.000390000000000;

		};
		~Config(){};
		Vector4d q0;
		
	};
}

#endif
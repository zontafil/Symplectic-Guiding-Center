//EM field interface. Must implement vector potential and magnetic field

#ifndef EMFIELD_H
#define EMFIELD_H

namespace EMFields{
	class EMField
	{
		public:
			EMField(Config::Config* config){};
			~EMField(){};

			virtual Vector3d A(Vector3d x) = 0;
			virtual Vector3d B(Vector3d x) = 0;
	
	};
}

#endif
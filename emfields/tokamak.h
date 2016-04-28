#ifndef TOKAMAK_H
#define TOKAMAK_H

#include "guidingfield.h"

namespace GuidingFields{
	class Tokamak : public GuidingFieldConfiguration{
		public:
			Tokamak(Config::Config* config): GuidingFieldConfiguration(config){
				B0 = config->B0;
				R0 = config->R0;
				kt = config->kt;
				q = config->q;
			};
			~Tokamak(){};

			Vector3d guiding_A(Vector3d x);
			Vector3d guiding_B(Vector3d x);

			double B0;
			double R0;
			double kt;
			double q;
	};	

	Vector3d Tokamak::guiding_A(Vector3d x){
		Vector3d ret;
		ret(0) = 0;
		ret(1) = -B0*R0*log((R0+x(0))/R0)*kt;
		ret(2) = B0/(2.*q*(R0+x(0)))*( 2.*R0*(R0+x(0))*log((R0+x(0))/R0) - 2.*R0*x(0) -2.*x(0)*x(0) - x(1)*x(1) );
		return ret;

	}
	Vector3d Tokamak::guiding_B(Vector3d x){
		Vector3d ret;
		ret(0) = -B0*x(1)/(q*(R0+x(0)));
		ret(1) = -((B0 *( R0 *x(0) *(-2.) -2.*x(0)*x(0) + x(1)*x(1)))/(2. * q *(R0 + x(0))*(R0 + x(0))));
		ret(2) = -B0*R0/(R0+x(0))*kt;
		return ret;

	}


}

#endif
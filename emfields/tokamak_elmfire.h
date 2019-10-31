//tokamak field configuration. (case B)
// add ELMFIRE poloidal component
// WARNING: this em field does not have the A potential, use only with integrators that don't need it

#ifndef TOKAMAK_ELMFIRE_H
#define TOKAMAK_ELMFIRE_H

#include "emField.h"

namespace EMFields{
	class TokamakElmfire : public EMField{
		private:
			double B0;
			double R0;
			double kt;
			double q;
			double aj;
			double a;
		public:
			TokamakElmfire(Config::Config* config): EMField(config), B0(config->B0), R0(config->R0), kt(config->kt), q(config->q) {
				this->aj = 1.5;
				this->a = 1.;
			};
			~TokamakElmfire(){};

			Vector3d A(Vector3d x);
			Vector3d B(Vector3d x);
	};	

	Vector3d TokamakElmfire::A(Vector3d x){
		Vector3d ret(0,0,0);
		return ret;
	}
	Vector3d TokamakElmfire::B(Vector3d x){
		Vector3d ret;

		double r = sqrt(x(0)*x(0) + x(1)*x(1));		
		double alpha = B0 / q / (R0 + x(0)) * (1./r/r * (1. - pow(1 - r*r/a/a, aj + 1)));

		ret(0) = - alpha * x(1);
		ret(1) = alpha * x(0);
		ret(1) += -(B0 * x(1)*x(1))/(2. * q *(R0 + x(0))*(R0 + x(0)));
		ret(2) = -B0*R0/(R0+x(0))*kt;
		return ret;

	}


}

#endif
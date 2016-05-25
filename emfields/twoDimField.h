//2D field configuration. (case A)
// see  6.2.2 paragraph

#ifndef TWODIMFIELD_H
#define TWODIMFIELD_H

#include "emField.h"

namespace EMFields{
	class TwoDimField : public EMField{
		public:
			TwoDimField(Config::Config* config): EMField(config){
				B0 = config->B0;
				R0 = config->R0;
				kt = config->kt;
				q = config->q;
			};
			~TwoDimField(){};

			Vector3d A(Vector3d x);
			Vector3d B(Vector3d x);

			double B0;
			double R0;
			double kt;
			double q;
	};	

	Vector3d TwoDimField::A(Vector3d x){
		Vector3d ret(0.,0.,0.);
		ret(0)=(-1.)*0.05/3. * x(1)*x(1)*x(1);
		ret(1)=0.05/3. * x(0)*x(0)*x(0)/4. +x(0);
		ret(2)=0;
		return ret;
	}
	Vector3d TwoDimField::B(Vector3d x){
		Vector3d ret(0.,0.,0.);
		ret(2)=1.+0.05*(x(0)*x(0)/4. + x(1)*x(1));
		return ret;
	}


}

#endif
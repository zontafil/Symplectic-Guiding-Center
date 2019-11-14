//force-free (perturbed and unperturbed) field configuration. (case C)
// see  6.2.3 paragraph

#ifndef FORCEFREE_H
#define FORCEFREE_H

#include <gsl/gsl_sf_bessel.h>
#include "emField.h"

#define BesselJ gsl_sf_bessel_Jn

namespace EMFields{
	class ForceFree : public EMField{
		private:
			Vector3d calc_B_nm(Vector3d x,int m,int n);
			Vector3d calc_A_nm(Vector3d x,int n,int m);
		public:
			ForceFree(Config::Config* config): EMField(config){
				B0 = config->B0;
				R0 = config->R0;
				kt = config->kt;
				q = config->q;

				forcefree_pert = config->forcefree_pert;
				forcefree_kff = config->forcefree_kff;
				a = config->forcefree_a;

			};
			~ForceFree(){};

			Vector3d A(Vector3d x);
			Vector3d B(Vector3d x);

			double B0;
			double R0;
			double kt;
			double q;

			bool forcefree_pert;
			double forcefree_kff, a;

	};	

	Vector3d ForceFree::calc_B_nm(Vector3d x,int m,int n){
	  return Vector3d
	  ((1/(a*(x(0)*x(0) + x(1)*x(1))))*((-m)*(a - n)*gsl_sf_bessel_Jn(m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*(x(1)*cos(n*x(2) + m*atan2(x(1),x(0))) + x(0)*sin(n*x(2) + m*atan2(x(1),x(0)))) + 
	    sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1))*gsl_sf_bessel_Jn(-1 + m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*(a*x(1)*cos(n*x(2) + m*atan2(x(1),x(0))) - n*x(0)*sin(n*x(2) + m*atan2(x(1),x(0))))), 
	  (1/(a*(x(0)*x(0) + x(1)*x(1))))*(m*(a - n)*gsl_sf_bessel_Jn(m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*(x(0)*cos(n*x(2) + m*atan2(x(1),x(0))) - x(1)*sin(n*x(2) + m*atan2(x(1),x(0)))) - 
	    sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1))*gsl_sf_bessel_Jn(-1 + m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*(a*x(0)*cos(n*x(2) + m*atan2(x(1),x(0))) + n*x(1)*sin(n*x(2) + m*atan2(x(1),x(0))))), 
	  ((a - n)*(a + n)*gsl_sf_bessel_Jn(m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*cos(n*x(2) + m*atan2(x(1),x(0))))/a);
	}

	Vector3d ForceFree::calc_A_nm(Vector3d x,int n,int m){
	  return Vector3d(((sqrt((a - n)*(a + n))*x(1)*gsl_sf_bessel_Jn(-1 + m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*cos(n*x(2) + m*atan2(x(1),x(0))))/sqrt(x(0)*x(0) + x(1)*x(1)) - 
	    (m*gsl_sf_bessel_Jn(m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*(x(1)*cos(n*x(2) + m*atan2(x(1),x(0))) + x(0)*sin(n*x(2) + m*atan2(x(1),x(0)))))/(x(0)*x(0) + x(1)*x(1)))/a, 
	  (-((sqrt((a - n)*(a + n))*x(0)*gsl_sf_bessel_Jn(-1 + m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*cos(n*x(2) + m*atan2(x(1),x(0))))/sqrt(x(0)*x(0) + x(1)*x(1))) + 
	    (m*gsl_sf_bessel_Jn(m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*(x(0)*cos(n*x(2) + m*atan2(x(1),x(0))) - x(1)*sin(n*x(2) + m*atan2(x(1),x(0)))))/(x(0)*x(0) + x(1)*x(1)))/a, 
	  gsl_sf_bessel_Jn(m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*cos(n*x(2) + m*atan2(x(1),x(0))));
	}

	Vector3d ForceFree::A(Vector3d x){
		if (forcefree_pert) return 1.*calc_A_nm(x,0,0)+forcefree_kff*calc_A_nm(x,1,1)+forcefree_kff*calc_A_nm(x,1,2);
		else return calc_A_nm(x,0,0);
	}
	Vector3d ForceFree::B(Vector3d x){
		if (forcefree_pert) return 1.*calc_B_nm(x,0,0)+forcefree_kff*calc_B_nm(x,1,1)+forcefree_kff*calc_B_nm(x,1,2);
		else return calc_B_nm(x,0,0);
	}


}

#endif
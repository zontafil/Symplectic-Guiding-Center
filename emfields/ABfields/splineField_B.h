//tokamak field configuration. (case B)
// add ELMFIRE poloidal component
// WARNING: this em field does not have the A potential, use only with integrators that don't need it

#ifndef SPLINEFIELDB_H
#define SPLINEFIELDB_H

#include "emField.h"
#include <stdexcept>
#include "../eqdskReader/eqdskReader.h"
#include "../eqdskReader/eqdskPsiInterp.h"
#include "../ascot5-spline/common_spline.h"
#include "../utils/coordinatesUtils.h"

namespace EMFields{
	class SplineField_B : public EMField{
		private:
            eqdsk eqdskObj;
            interp2D_data* psi_spline_c;
            interp1D_data* fpol_spline_c;
		public:
			SplineField_B(Config::Config* config);
			~SplineField_B(){};

			Vector3d A(Vector3d x);
			Vector3d B(Vector3d x);
	};	

    SplineField_B::SplineField_B(Config::Config* config): EMField(config) {
        if (config->eqdsk_file.empty()) {
            throw invalid_argument("EQDSK input file unset. Set config->eqdsk_file");
        }
        eqdskObj = readEqdskFile(config->eqdsk_file);
        psi_spline_c = eqdskPsiInterp(eqdskObj);
        fpol_spline_c = eqdskFpolInterp(eqdskObj);

        cout << "EQDSK: range r: " << eqdskObj.r_min << " " << eqdskObj.r_max << endl;
        cout << "EQDSK: range z: " << eqdskObj.z_min << " " << eqdskObj.z_max << endl;
	};

	Vector3d SplineField_B::A(Vector3d x){
		Vector3d ret(0,0,0);
		return ret;
	}
	Vector3d SplineField_B::B(Vector3d x){
        realnum r = sqrt(x[0]*x[0] + x[1]*x[1]);

        // compute B, dB in cyl coordinates
        BdB_rz BdB_cyl = evalBrz(r, x[2], psi_spline_c, fpol_spline_c, eqdskObj);
        Vector3d Bcyl(BdB_cyl.BR, BdB_cyl.Bp, BdB_cyl.Bz);

        // build B, B and |B| (cartesian)
        Vector3d B = cyl2cart(Bcyl, x);

        return B;
	}


}

#endif
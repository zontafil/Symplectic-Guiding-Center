#ifndef SPLINEFIELD_H
#define SPLINEFIELD_H

#include "guidingfield.h"
#include <stdexcept>
#include "./eqdskReader/eqdskReader.h"
#include "./eqdskReader/eqdskPsiInterp.h"
#include "ascot5-spline/common_spline.h"

namespace EMFields{

    Vector3d cyl2cart(Vector3d v, Vector3d x) {
        realnum r = sqrt(x[0]*x[0] + x[1]*x[1]);
        Vector3d ret;
        ret[0] = (v[0]*x[0] - v[1]*x[1]) / r;
        ret[1] = (v[0]*x[1] + v[1]*x[0]) / r;
        ret[2] = v[2];
        return ret;
    }

	template <int DIM> class SplineField: public GuidingFieldConfiguration<DIM>
	{
        private:
            eqdsk eqdskObj;
		public:
			SplineField(Config::Config* config);
			~SplineField(){};

			//compute the field from q
			GuidingField compute(Matrix<double,DIM/2,1> q);
	};

	template <int DIM> SplineField<DIM>::SplineField(Config::Config* config): GuidingFieldConfiguration<DIM>(config) {
        if (config->eqdsk_file.empty()) {
            throw invalid_argument("EQDSK input file unset. Set config->eqdsk_file");
        }
        eqdskObj = readEqdskFile(config->eqdsk_file);

        cout << "EQDSK: range r: " << eqdskObj.r_min << " " << eqdskObj.r_max << endl;
        cout << "EQDSK: range z: " << eqdskObj.z_min << " " << eqdskObj.z_max << endl;
	};

    template <int DIM> GuidingField SplineField<DIM>::compute(Matrix<double,DIM/2,1> q){
        GuidingField ret;

        interp2D_data* psi_spline_c = eqdskPsiInterp(eqdskObj);
        interp1D_data* fpol_spline_c = eqdskFpolInterp(eqdskObj);

        // TODO
        Vector3d x = q.head(3);
        realnum r = sqrt(x[0]*x[0] + x[1]*x[1]);

        // compute B, dB in cyl coordinates
        BdB_rz BdB_cyl = evalBrz(r, x[2], psi_spline_c, fpol_spline_c, eqdskObj);
        Vector3d Bcyl(BdB_cyl.BR, BdB_cyl.Bp, BdB_cyl.Bz);

        ret.B = cyl2cart(Bcyl, x);
        ret.b = ret.B.normalized();

        // build curl B and Bdag
        Vector3d Bcurl_cyl;
        Bcurl_cyl[0] = BdB_cyl.dBz_dp / r - BdB_cyl.dBp_dz;
        Bcurl_cyl[1] = BdB_cyl.dBR_dz - BdB_cyl.dBz_dR;
        Bcurl_cyl[2] = Bcyl[1] / r + BdB_cyl.dBp_dR - BdB_cyl.dBp_dp / r;
        Vector3d Bcurl = cyl2cart(Bcurl_cyl, x);
        ret.Bdag = ret.B + q[3] * Bcurl;

        // build grad|B|
        realnum Bnorm = ret.B.norm();
        Vector3d gradB_cyl;
        gradB_cyl[0] = Bcyl.dot(Vector3d(BdB_cyl.dBR_dR, BdB_cyl.dBp_dR, BdB_cyl.dBz_dR));
        gradB_cyl[1] = Bcyl.dot(Vector3d(BdB_cyl.dBR_dp, BdB_cyl.dBp_dp, BdB_cyl.dBz_dp));
        gradB_cyl[2] = Bcyl.dot(Vector3d(BdB_cyl.dBR_dz, BdB_cyl.dBp_dz, BdB_cyl.dBz_dz));
        gradB_cyl /= Bnorm;
        ret.B_grad = cyl2cart(gradB_cyl, x);

        return ret;
    }


}

#endif
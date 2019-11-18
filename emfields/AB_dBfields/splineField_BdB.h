#ifndef SPLINEFIELD_H
#define SPLINEFIELD_H

#include "../AB_dB_Field.h"
#include <stdexcept>
#include "../eqdskReader/eqdskReader.h"
#include "../eqdskReader/eqdskPsiInterp.h"
#include "../ascot5-spline/common_spline.h"
#include "../utils/coordinatesUtils.h"

namespace EMFields{

	template <int DIM> class SplineField_BdB: public AB_dB_FieldBuilder<DIM>
	{
        private:
            eqdsk eqdskObj;
            interp2D_data* psi_spline_c;
            Spline1D fpol_spline;
		public:
			SplineField_BdB(Config::Config* config);
			~SplineField_BdB(){
                free(psi_spline_c->c);
                free(psi_spline_c);
            };

			//compute the field from q
			GuidingField compute(Matrix<double,DIM/2,1> q);
	};

	template <int DIM> SplineField_BdB<DIM>::SplineField_BdB(Config::Config* config): AB_dB_FieldBuilder<DIM>(config) {
        if (config->eqdsk_file.empty()) {
            throw invalid_argument("EQDSK input file unset. Set config->eqdsk_file");
        }
        eqdskObj = readEqdskFile(config->eqdsk_file);
        psi_spline_c = eqdskPsiInterp(eqdskObj);
        fpol_spline = eqdskFpolInterpEigen(eqdskObj);

        cout << "EQDSK: range r: " << eqdskObj.r_min << " " << eqdskObj.r_max << endl;
        cout << "EQDSK: range z: " << eqdskObj.z_min << " " << eqdskObj.z_max << endl;
        cout << "EQDSK: range psi: " << eqdskObj.simag << " " << eqdskObj.sibry << endl;

        // DELME: compute B field around a point
        // Vector3d x0(0.31266191772552109907, 1.40681096836614916157, 0.26799155264395102538);
        // for (double r = 1.2; r <= 2.; r+=0.001) {
        //     Vector4d q(r,0,x0[2],1.);

        //     // realnum r = sqrt(x*x + x0[1]*x0[1]);
        //     BdB_rz BdB_cyl = evalBrz(r, x0[2], psi_spline_c, fpol_spline, eqdskObj);

        //     // compute B, dB in cyl coordinates
        //     GuidingField BdB = this->compute(q);
        //     cout << r << " ";
        //     cout << BdB_cyl.fpol_comp1 << " "; //fpol_i_r
        //     cout << BdB_cyl.psi_ir << " ";
        //     cout << BdB_cyl.psi_dr << " "; // 4
        //     cout << BdB_cyl.fpol << " ";
        //     cout << BdB_cyl.psi << " ";
        //     cout << BdB_cyl.psi_comp1 << " ";
        //     cout << BdB_cyl.dpsi_dR << " "; //8
        //     cout << BdB_cyl.dpsi_dz << " ";
        //     cout << BdB_cyl.d2psi_dR2 << " ";
        //     cout << BdB_cyl.d2psi_dz2 << " ";
        //     cout << BdB_cyl.d2psi_dRdz << " ";
        //     cout << BdB_cyl.dfpol_dpsi << " "; // 13
        //     cout << BdB_cyl.d2fpol_dpsi2 << " "; // 14
        //     cout << BdB_cyl.dpsi_dR << " ";
        //     cout << BdB_cyl.BR << " ";
        //     cout << BdB_cyl.Bp << " ";
        //     cout << BdB_cyl.Bz << " ";
        //     cout << BdB_cyl.dBR_dR << " ";
        //     cout << BdB_cyl.dBR_dp << " "; // 20
        //     cout << BdB_cyl.dBR_dz << " ";
        //     cout << BdB_cyl.dBp_dR << " ";
        //     cout << BdB_cyl.dBp_dp << " ";
        //     cout << BdB_cyl.dBp_dz << " ";
        //     cout << BdB_cyl.dBz_dR << " ";
        //     cout << BdB_cyl.dBz_dp << " ";
        //     cout << BdB_cyl.dBz_dz << " ";
        //     cout << BdB.B.transpose() << " "; //28
        //     cout << BdB.B_grad.transpose() << " "; //31
        //     cout << BdB.GradB_cyl.transpose() << " "; //34
        //     cout << BdB.Bcyl.transpose() << " "; //37
        //     cout << "0 "; //40
        //     cout << BdB.Bnorm << " "; // 41
        //     cout << BdB.Bdag.transpose() << " "; // 42
        //     cout << endl;
        // }


	};

    template <int DIM> GuidingField SplineField_BdB<DIM>::compute(Matrix<double,DIM/2,1> q){
        GuidingField ret;

        Vector3d x = q.head(3);
        realnum r = sqrt(x[0]*x[0] + x[1]*x[1]);

        // compute B, dB in cyl coordinates
        BdB_rz BdB_cyl = evalBrz(r, x[2], psi_spline_c, fpol_spline, eqdskObj);
        Vector3d Bcyl(BdB_cyl.BR, BdB_cyl.Bp, BdB_cyl.Bz);

        // build B, B and |B| (cartesian)
        ret.B = cyl2cart(Bcyl, x);
        ret.b = ret.B.normalized();
        ret.Bnorm = ret.B.norm();

        // build curl B and Bdag (cyl and cartesian)
        Vector3d Bcurl_cyl;
        Bcurl_cyl[0] = BdB_cyl.dBz_dp / r - BdB_cyl.dBp_dz;
        Bcurl_cyl[1] = BdB_cyl.dBR_dz - BdB_cyl.dBz_dR;
        Bcurl_cyl[2] = Bcyl[1] / r + BdB_cyl.dBp_dR - BdB_cyl.dBp_dp / r;
        Vector3d Bcurl = cyl2cart(Bcurl_cyl, x);
        ret.Bdag = ret.B + q[3] * Bcurl / ret.Bnorm;

        // build grad|B| (cyl and cartesian)
        Vector3d gradB_cyl;
        gradB_cyl[0] = Bcyl.dot(Vector3d(BdB_cyl.dBR_dR, BdB_cyl.dBp_dR, BdB_cyl.dBz_dR));
        gradB_cyl[1] = Bcyl.dot(Vector3d(BdB_cyl.dBR_dp, BdB_cyl.dBp_dp, BdB_cyl.dBz_dp)) / r;
        gradB_cyl[2] = Bcyl.dot(Vector3d(BdB_cyl.dBR_dz, BdB_cyl.dBp_dz, BdB_cyl.dBz_dz));
        gradB_cyl /= ret.Bnorm;
        ret.GradB_cyl = gradB_cyl;
        ret.Bcyl = Bcyl;
        ret.B_grad = cyl2cart(gradB_cyl, x);

        return ret;
    }


}

#endif
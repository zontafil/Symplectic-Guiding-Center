#include "./eqdskReader.h"
#include "../ascot5-spline/inter2Dexpl.h"
#include <algorithm> 
#include "./eqdskPsiInterp.h"
#include <stdexcept>
#include <iostream>

using namespace std;

interp2D_data* eqdskPsiInterp(eqdsk eqdsk_data) {
    interp2D_data* ret = new interp2D_data();
    realnum* f = eqdsk_data.psirz;
    int nr = eqdsk_data.nr;
    int nz = eqdsk_data.nz;
    realnum r_min = eqdsk_data.r_min;
    realnum r_max = eqdsk_data.r_max;
    realnum r_grid = eqdsk_data.r_grid;
    realnum z_min = eqdsk_data.z_min;
    realnum z_max = eqdsk_data.z_max;
    realnum z_grid = eqdsk_data.z_grid;

    a5err err;
    if ((err = interp2Dexpl_init(ret, f, nr, nz, r_min, r_max, r_grid, z_min, z_max, z_grid))) {
        throw invalid_argument("2D psi(r,z) fitting failed: "+ to_string(err));
    }

    return ret;
};
Spline1D eqdskFpolInterpEigen(eqdsk eqdsk_data) {
        int n_r = eqdsk_data.nr;
        Eigen::RowVectorXd fpol, psi;
        fpol.resize(n_r);
        psi.resize(eqdsk_data.nr);

        for (int i = 0; i < n_r; i++) {
            fpol(i) = eqdsk_data.fpol[i];
            psi(i) = i;
        }

        Spline1D spline = SplineFitting1D::Interpolate(fpol, 4, psi);
        return spline;
}

int eqdskFpolEvalEigen(realnum *fdf, eqdsk eqdsk_data, Spline1D spline, realnum psi) {
    realnum r_min = eqdsk_data.simag;
    realnum r_max = eqdsk_data.sibry;
    realnum psi_normalized = (psi - r_min) / (r_max - r_min) * double(eqdsk_data.nr);

    const auto eval = spline.derivatives(psi_normalized, 2);

    fdf[0] = eval[0];
    fdf[1] = eval[1] / (r_max - r_min) * eqdsk_data.nr;
    fdf[2] = eval[2] / pow(r_max - r_min, 2) * pow(eqdsk_data.nr, 2);
    fdf[3] = psi_normalized;

    return 0;
}

BdB_rz evalBrz(realnum r, realnum z, interp2D_data* interp2Dc, Spline1D splineFpol, eqdsk eqdskData) {
    a5err err;
    BdB_rz ret;

    // interpolate psi, dpsi
    realnum psi_dpsi[6];
    if ((err = interp2Dexpl_eval_dB(psi_dpsi, interp2Dc, r, z))) {
        throw invalid_argument("psi Eval error. Point outside interpolation region");
    }


    realnum psi = psi_dpsi[0];
    realnum dpsi_dR = psi_dpsi[1];
    realnum dpsi_dz = psi_dpsi[2];
    realnum d2psi_dR2 = psi_dpsi[3];
    realnum d2psi_dz2 = psi_dpsi[4];
    realnum d2psi_dRdz = psi_dpsi[5];

    // interpolate fpol, if inside the main plasma boundary
    realnum fpol;
    realnum dfpol_dpsi;
    realnum d2fpol_dpsi2;

    realnum fpol_df[3];
    if ((psi < max(eqdskData.sibry, eqdskData.simag)) && (psi > min(eqdskData.sibry, eqdskData.simag))
        && (z < eqdskData.sepmaxz) && (z > eqdskData.sepminz)) {
        // most likely in the main plasma

        // if ((err = interp1Dcomp_eval_dB(fpol_df, interp1Dc, psi))) {
        if ((err = eqdskFpolEvalEigen(fpol_df, eqdskData, splineFpol, psi))) {
            throw invalid_argument("fpol Eval error " + to_string(err));
        }
        fpol = fpol_df[0];
        dfpol_dpsi = fpol_df[1];
        d2fpol_dpsi2 = fpol_df[2];
    } else {
        cout << "WARNING! OUTSIDE MAIN PLASMA " << endl;
        // most likely outside the main plasma

        fpol = eqdskData.fpol[eqdskData.nr - 1];
        dfpol_dpsi = 0.;
        d2fpol_dpsi2 = 0.;
    }

    // evaluate the magnetic field
    ret.BR = -dpsi_dz / r;
    ret.Bp = fpol / r;
    ret.Bz = dpsi_dR / r;

    // evaluate the derivatives
    ret.dBR_dR = dpsi_dz/(r*r)-d2psi_dRdz/r;
    ret.dBR_dp = 0.;
    ret.dBR_dz = -d2psi_dz2/r;
    ret.dBp_dR = -fpol/(r*r)+dfpol_dpsi*dpsi_dR/r;
    ret.dBp_dp = 0.;
    ret.dBp_dz = dfpol_dpsi*dpsi_dz/r;
    ret.dBz_dR = -dpsi_dR/(r*r) + d2psi_dR2/r;
    ret.dBz_dp = 0.;
    ret.dBz_dz = d2psi_dRdz/r;
    ret.psi = psi;
    ret.dpsi_dR = dpsi_dR;
    ret.dpsi_dz = dpsi_dz;
    ret.d2psi_dR2 = d2psi_dR2;
    ret.d2psi_dz2 = d2psi_dz2;
    ret.d2psi_dRdz = d2psi_dRdz;

    ret.fpol = fpol;
    ret.dfpol_dpsi = dfpol_dpsi;
    ret.d2fpol_dpsi2 = d2fpol_dpsi2;
    ret.dpsi_dR = dpsi_dR;

    return ret;
}
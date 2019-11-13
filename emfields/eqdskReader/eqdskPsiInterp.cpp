#include "./eqdskReader.h"
#include "../ascot5-spline/inter2Dexpl.h"
#include "../ascot5-spline/interp1Dcomp.h"
#include <algorithm> 
#include "./eqdskPsiInterp.h"
#include <stdexcept>

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
        throw invalid_argument("Invalid Hamiltonian system "+ to_string(err));
    }

    return ret;
};
interp1D_data* eqdskFpolInterp(eqdsk eqdsk_data) {
    interp1D_data* ret = new interp1D_data();
    realnum r_min = eqdsk_data.simag;
    realnum r_max = eqdsk_data.sibry;
    int n_r = eqdsk_data.nr;
    realnum r_grid = (r_max - r_min) / (n_r - 1);
    realnum* f = eqdsk_data.fpol;

    int err;
    if ((err = interp1Dcomp_init(ret, f, n_r, r_min, r_max, r_grid))) {
        throw invalid_argument("Interpolation 1D error " + to_string(err));
    }

    return ret;
};

BdB_rz evalBrz(realnum r, realnum z, interp2D_data* interp2Dc, interp1D_data* interp1Dc, eqdsk eqdskData) {
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
    realnum d2psi_dRdz = psi_dpsi[5];

    // interpolate fpol, if inside the main plasma boundary
    realnum fpol;
    realnum dfpol_dpsi;

    if ((psi < max(eqdskData.sibry, eqdskData.simag)) && (psi > min(eqdskData.sibry, eqdskData.simag))
        && (z < eqdskData.sepmaxz) && (z > eqdskData.sepminz)) {
        // most likely in the main plasma

        realnum fpol_df[3];
        if ((err = interp1Dcomp_eval_dB(fpol_df, interp1Dc, psi))) {
            throw invalid_argument("fpol Eval error " + to_string(err));
        }
        fpol = fpol_df[0];
        dfpol_dpsi = fpol_df[1];
    } else {
        // most likely outside the main plasma

        fpol = eqdskData.fpol[eqdskData.nr - 1];
        dfpol_dpsi = 0.;
    }

    // evaluate the magnetic field
    ret.BR = -dpsi_dz / r;
    ret.Bp = fpol / r;
    ret.Bz = dpsi_dR / r;

    // evaluate the derivatives
    ret.dBR_dR = dpsi_dz/(r*r)-d2psi_dRdz/r;
    ret.dBR_dp = 0.;
    ret.dBR_dz = -d2psi_dRdz/r;
    ret.dBp_dR = -fpol/(r*r)+dfpol_dpsi*dpsi_dR/r;
    ret.dBp_dp = 0.;
    ret.dBp_dz = dfpol_dpsi*dpsi_dz/r;
    ret.dBz_dR = -dpsi_dR/(r*r) + d2psi_dR2/r;
    ret.dBz_dp = 0.;
    ret.dBz_dz = d2psi_dRdz/r;

    return ret;
}
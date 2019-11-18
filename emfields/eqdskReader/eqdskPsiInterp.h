#ifndef BFIELD_COMPUTE_H
#define BFIELD_COMPUTE_H

#include "../ascot5-spline/inter2Dexpl.h"
#include "./eqdskReader.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/Splines>

typedef Eigen::Spline<double, 1, 4> Spline1D;
typedef Eigen::SplineFitting<Spline1D> SplineFitting1D;

typedef struct {
    realnum BR;
    realnum Bp;
    realnum Bz;
    realnum dBR_dR;
    realnum dBR_dp;
    realnum dBR_dz;
    realnum dBp_dR;
    realnum dBp_dp;
    realnum dBp_dz;
    realnum dBz_dR;
    realnum dBz_dp;
    realnum dBz_dz;
    realnum fpol;
    realnum dfpol_dpsi;
    realnum d2fpol_dpsi2;
    realnum psi;
    realnum dpsi_dR;
    realnum dpsi_dz;
    realnum d2psi_dR2;
    realnum d2psi_dz2;
    realnum d2psi_dRdz;
} BdB_rz;

Spline1D eqdskFpolInterpEigen(eqdsk eqdsk_data);
interp2D_data* eqdskPsiInterp(eqdsk eqdsk_data);
BdB_rz evalBrz(realnum r, realnum z, interp2D_data* interp2Dc, Spline1D splineFpol, eqdsk eqdskData);


#endif
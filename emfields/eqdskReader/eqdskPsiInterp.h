#ifndef BFIELD_COMPUTE_H
#define BFIELD_COMPUTE_H

#include "../ascot5-spline/inter2Dexpl.h"
#include "../ascot5-spline/interp1Dcomp.h"
#include "./eqdskReader.h"

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
} BdB_rz;

interp2D_data* eqdskPsiInterp(eqdsk eqdsk_data);
interp1D_data* eqdskFpolInterp(eqdsk eqdsk_data);
BdB_rz evalBrz(realnum r, realnum z, interp2D_data* interp2Dc, interp1D_data* interp1Dc, eqdsk eqdskData);

#endif
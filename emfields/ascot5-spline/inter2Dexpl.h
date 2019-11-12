#ifndef INTERP2DEXP_H
#define INTERP2DEXP_H

#include "./common_spline.h"

int interp2Dexpl_init(interp2D_data* str, realnum* f, int n_r, int n_z,
                       realnum r_min, realnum r_max, realnum r_grid,
                       realnum z_min, realnum z_max, realnum z_grid);
a5err interp2Dexpl_eval_dB(realnum* B_dB, interp2D_data* str, realnum r, realnum z);
a5err interp2Dexpl_eval_dB(realnum* B_dB, interp2D_data* str, realnum r, realnum z);



#endif
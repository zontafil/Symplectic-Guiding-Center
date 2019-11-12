#ifndef INTERP1DCOMP_H
#define INTERP1DCOMP_H

#include "./common_spline.h"

int interp1Dcomp_init(interp1D_data* str, realnum* f, int n_r,
                      realnum r_min, realnum r_max, realnum r_grid);
int interp1Dcomp_eval_B(realnum* B, interp1D_data* str, realnum r);
int interp1Dcomp_eval_dB(realnum* B_dB, interp1D_data* str, realnum r);


#endif
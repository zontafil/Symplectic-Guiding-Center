/**
 * @file spline.h
 * @brief Header file for spline_expl.c and spline_comp.c
 */
#ifndef SPLINE_H
#define SPLINE_H
#include "common_spline.h"

void spline(realnum* f, int n, int bc, realnum* c);
void splinecomp(realnum* f, int n, int bc, realnum* c);
void spline1Dcomp(realnum* f, int n, int bc, realnum* c);
#endif

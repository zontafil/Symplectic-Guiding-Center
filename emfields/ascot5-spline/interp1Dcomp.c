/**
 * @file interp1Dcomp.c
 * @brief Bicubic spline interpolation in compact form
 */
#include <math.h>
#include "./interp1Dcomp.h"
#include <stdlib.h>
#include "./spline.h"

/**
 * @brief Calculate bicubic spline interpolation coefficients for scalar 1D data
 *
 * This function calculates the bicubic spline interpolation coefficients for
 * the given data and stores them in the data struct. Compact  cofficients are
 * calculated directly.
 *
 * @todo Error checking
 *
 * @param str data struct for data interpolation
 * @param f 1D data to be interpolated
 * @param n_r number of data points in the r direction
 * @param r_min minimum value of the r axis
 * @param r_max maximum value of the r axis
 */
int interp1Dcomp_init(interp1D_data* str, realnum* f, int n_r,
                      realnum r_min, realnum r_max, realnum r_grid) {
    int err = 0;

    /* Initialize and fill the data struct */
    str->n_r = n_r;
    str->r_min = r_min;
    str->r_max = r_max;
    str->r_grid = r_grid;
    str->c = (realnum*)malloc(n_r*2*sizeof(realnum));

    if(str->c == NULL) {
        err = 1;
    }
    else {
        /* Note! Explicitly using non-periodic boundary conditions. */
        spline1Dcomp(f,str->n_r,0,str->c);
    }

    return err;
}

/**
 * @brief Evaluate interpolated value of 1D scalar field
 *
 * This function evaluates the interpolated value of a 1D scalar field using
 * bicubic spline interpolation coefficients of the compact form.
 *
 * @todo Check discrepency to ascot4 and explicit version
 * @todo Error checking
 *
 * @param B variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param r r-coordinate
 */
int interp1Dcomp_eval_B(realnum* B, interp1D_data* str, realnum r) {
    int i_r = (r-str->r_min)/str->r_grid; /**< index for r variable */
    realnum dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
                                                               current cell */
    realnum dri = 1.0-dr;
    realnum dri3 = dri*(dri*dri-1.0);
    realnum rg2 = str->r_grid*str->r_grid;        /**< Square of cell length in r direction */

    int n = i_r*2;       /**< Index jump to cell */
    int r1 = 2;                         /**< Index jump one r forward */

    int err = 0;

    /* Check that the point is not outside the evaluation regime */
    if(r < str->r_min || r > str->r_max) {
        err = 1;
    }
    else {
        *B = (
              dri*str->c[n] +
              dr*str->c[n+r1] +
              (rg2/6)*(dri3*str->c[n+1] + dri3*str->c[n+r1+1])
              );
    }
    return err;
}

/**
 * @brief Evaluate interpolated value of 1D scalar field and its 1st and 2nd derivatives
 *
 * This function evaluates the interpolated value of a 1D scalar field and
 * its 1st and 2nd derivatives using bicubic spline interpolation coefficients
 * of the compact form.
 *
 * @todo Check discrepency to ascot4 and explicit version
 * @todo Error checking
 *
 * @param B_dB array in which to place the evaluated values
 * @param str data struct for data interpolation
 * @param r r-coordinate
 */
int interp1Dcomp_eval_dB(realnum* B_dB, interp1D_data* str, realnum r) {
    int i_r = (r-str->r_min)/str->r_grid;                   /**< index for r variable */
    realnum dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
                                                               current cell */
    realnum dr3dr = 3*dr*dr-1;           /**< r-derivative of dr3, not including 1/r_grid */
    realnum dri = 1.0-dr;
    realnum dri3 = dri*(dri*dri-1);
    realnum dri3dr = -3*dri*dri+1;       /**< r-derivative of dri3, not including 1/r_grid */
    realnum rg = str->r_grid;            /**< Cell length in r direction */
    realnum rg2 = rg*rg;
    realnum rgi = 1.0/rg;

    int n = i_r*2;     /**< Index jump to cell */
    int r1 = 2;                       /**< Index jump one r forward */

    int err = 0;

    /* Check that the point is not outside the evaluation regime */
    if(r < str->r_min || r > str->r_max) {
        err = 1;
    }
    else {
        /* f */
        B_dB[0] = (
                   dri*str->c[n] +
                   dr*str->c[n+r1] +
                   (rg2/6)*(dri3*str->c[n+1] + dri3*str->c[n+r1+1])
                   );

        /* df/dr */
        B_dB[1] = (rgi*(str->c[n+r1] - str->c[n])  +
                   (rg/6)*(dr3dr*str->c[n+r1+1] +dri3dr*str->c[n+1])
                   );

        /* d2f/dr2 */
        B_dB[2] = (
                   dri*str->c[n+1] + dr*str->c[n+r1+1]
                   );
    }
    return err;
}
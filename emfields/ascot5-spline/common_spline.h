/**
 * @file common_spline.h
 * @brief Header file for common stuff (not included in ascot5)
 */
#ifndef COMMON_SPLINE_H
#define COMMON_SPLINE_H

typedef double realnum;   /**< Singe precision float    */
typedef unsigned long int a5err;

/**
* @brief Boundary conditions for the spline interpolation.
*/
enum boundaryCondition {
    NATURALBC  = 0, /**< Second derivative is zero                            */
    PERIODICBC = 1  /**< Function has same value and derivatives on both ends */
};

typedef struct {
    int n_r;                  /**< number of r grid points */
    realnum r_min;               /**< minimum r coordinate in the grid */
    realnum r_max;               /**< r grid interval (r_max-r_min)/(n_r-1) */
    realnum r_grid;              /**< r grid interval (r_max-r_min)/(n_r-1) */
    realnum* c;                  /**< pointer to array with spline coefficients */
} interp1D_data;

typedef struct {
    int n_r;     /**< number of x grid points                        */
    int n_z;     /**< number of y grid points                        */
    int bc_x;    /**< boundary condition for x coordinate            */
    int bc_y;    /**< boundary condition for y coordinate            */
    realnum r_min;  /**< minimum x coordinate in the grid               */
    realnum r_max;  /**< maximum x coordinate in the grid               */
    realnum r_grid; /**< interval between two adjacent points in x grid */
    realnum z_min;  /**< minimum y coordinate in the grid               */
    realnum z_max;  /**< maximum y coordinate in the grid               */
    realnum z_grid; /**< interval between two adjacent points in y grid */
    realnum* c;     /**< pointer to array with spline coefficients      */
} interp2D_data;


#endif

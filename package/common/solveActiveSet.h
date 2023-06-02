/**
 * Copyright (C) Till Blaha 2022-2023
 * MAVLab -- Faculty of Aerospace Engineering -- Delft University of Techology
 */

/**
 * @file solveActiveSet.h
 * 
 * @brief Declares variants of the active-set solver, gateway funcs and utility
*/

#ifndef SOLVEACTIVESET_H
#define SOLVEACTIVESET_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __FAST_MATH__
#error "Do not use -ffast-math, needed for protection against nan"
#endif

#ifndef AS_N_U
#warning "AS_N_U not defined. Assuming 20"
#define AS_N_U 20
#endif

#ifndef AS_N_V
#warning "AS_N_V not defined. Assuming 6"
#define AS_N_V 6
#endif

#define AS_N_C (AS_N_U+AS_N_V)

#if defined(RECORD_COST) && !defined RECORD_COST_N
#define RECORD_COST_N 15
#endif

#if defined(RECORD_COST) && !defined RTOL
#define RTOL 1e-7
#endif

#if defined(RECORD_COST) && !defined CTOL
#define CTOL 1e-7
#endif

#ifdef USE_SINGLE_FLOAT
typedef float num_t;
#define TOL 1e-4
#else
typedef double num_t;
#define TOL 1e-7
#endif

typedef enum {
  ALLOC_SUCCESS = 0,
  ALLOC_ITER_LIMIT = 1,
  ALLOC_COST_BELOW_TOL = 2,
  ALLOC_COST_PLATEAU = 3,
  ALLOC_NAN_FOUND_Q = 4,
  ALLOC_NAN_FOUND_US = 5,
  } activeSetExitCode;

typedef enum {
  QR_NAIVE = 0, // mostly previous PPRZ
  QR = 1,
  CHOL = 2,
  CG = 3,
  } activeSetAlgoChoice;

/**
 * @brief Function type definition for any active-set solver
 * 
 * Solve min( ||Au-b||_2 )
 *          s.t. umin[i] <= u[i] <= umax[i]
 * 
 * @param A_col col-major representation of A
 * @param b rhs of equation
 * @param umin lower bound on u
 * @param umax upper bound on u
 * @param us initial guess. Becomes solution on exit, unless return > 0
 * @param Ws inital working set. Becomes working-set of iterate us on exit.
 * @param imax max number of iterations, default 100 if imax = 0
 * @param n_u Length of us
 * @param n_v Such that n_u+n_v is length of b
 * @param iter On exit: number of active-set iterations performed
 * @param n_free On exit: number of free variables at iterate us
 * @param costs if RECORD_COST defined: cost at first RECORD_COST_N iterations
 * 
 * @return 0 if success us is true solution, >0 if error. See activeSetExitCode
 */
typedef activeSetExitCode (*activeSetAlgo)(
  const num_t A_col[AS_N_C*AS_N_U], const num_t b[AS_N_C],
  const num_t umin[AS_N_U], const num_t umax[AS_N_U], num_t us[AS_N_U],
  int8_t Ws[AS_N_U], unsigned int imax, const int n_u, const int n_v,
  int *iter, int *n_free, num_t costs[RECORD_COST_N]);

/**
 * @brief Gateway function the allows switching 
 * 
 * Usage example:
 * activeSetExitCode alloc_result;
 * alloc_result = solveActiveSet(QR)(A_col, b, ...);
 * 
 * @param choice Which algorithm return
 * 
 * @return pointer to a activeSetAlgo
 */
activeSetAlgo solveActiveSet(activeSetAlgoChoice choice) {
    switch (choice) {
        case QR_NAIVE:
            return &solveActiveSet_qr_naive;
        case CHOL:
            return &solveActiveSet_chol;
        case CG:
            return &solveActiveSet_cg;
        default:
            return &solveActiveSet_qr;
    }
};

activeSetExitCode solveActiveSet_qr_naive(
  const num_t A_col[AS_N_C*AS_N_U], const num_t b[AS_N_C],
  const num_t umin[AS_N_U], const num_t umax[AS_N_U], num_t us[AS_N_U],
  int8_t Ws[AS_N_U], unsigned int imax, const int n_u, const int n_v,
  int *iter, int *n_free, num_t costs[]);
activeSetExitCode solveActiveSet_qr(
  const num_t A_col[AS_N_C*AS_N_U], const num_t b[AS_N_C],
  const num_t umin[AS_N_U], const num_t umax[AS_N_U], num_t us[AS_N_U],
  int8_t Ws[AS_N_U], unsigned int imax, const int n_u, const int n_v,
  int *iter, int *n_free, num_t costs[]);
activeSetExitCode solveActiveSet_chol(
  const num_t A_col[AS_N_C*AS_N_U], const num_t b[AS_N_C],
  const num_t umin[AS_N_U], const num_t umax[AS_N_U], num_t us[AS_N_U],
  int8_t Ws[AS_N_U], unsigned int imax, const int n_u, const int n_v,
  int *iter, int *n_free, num_t costs[]);
activeSetExitCode solveActiveSet_cg(
  const num_t A_col[AS_N_C*AS_N_U], const num_t b[AS_N_C],
  const num_t umin[AS_N_U], const num_t umax[AS_N_U], num_t us[AS_N_U],
  int8_t Ws[AS_N_U], unsigned int imax, const int n_u, const int n_v,
  int *iter, int *n_free, num_t costs[]);

#ifdef RECORD_COST
/**
 * @brief Compute penalty function cost as ||Au - b||^2
 * 
 * Uses sparse multiplication to save time.
 * 
 * @param A_col col-major representation of A
 * @param b rhs of equation
 * @param u test point u
 * @param n_u number of elements in u
 * @param n_v Such that n_u + n_v is length of b
 * 
 * @return cost of quadratic penalty function at point u
 */
num_t calc_cost(const num_t A_col[AS_N_C*AS_N_U], const num_t b[AS_N_C],
  const num_t u[AS_N_U], const int n_u, const int n_v)
{
	// checking cost in n_v*n_u+n_u time
	num_t cost = 0.;
	for (int i=0; i<(n_u+n_v); i++) {
		num_t i_cost = -b[i];
		if (i < n_v) {
			// dense part
			for (int j=0; j<n_u; j++)
				i_cost += A_col[i + j*(n_u+n_v)]*u[j];
		} else {
			// sparse part
			i_cost += A_col[i + (i-n_v)*(n_u+n_v)]*u[i-n_v]; 
		}
		cost += i_cost*i_cost;
	}
	return cost;
};
#endif

#endif


/* Mathilda — Taylor-series coefficients of an ODE solution (MPFR).
 *
 * For a first-order system  Y' = F(t, Y)  whose right-hand side is built from
 * + , * , powers with a constant exponent, Exp, Log, Sin, Cos, Tan, Sinh, Cosh,
 * Tanh and ArcTan, the Taylor coefficients of the solution through (t0, Y0)
 * are computed exactly (up to rounding) by automatic differentiation in Taylor
 * mode: the RHS expression trees are flattened once into a tape of elementary
 * operations, and the k-th coefficient of every tape node is obtained from the
 * lower-order ones by the standard recurrences (Cauchy product, the
 * a^alpha / exp / log / sin-cos recurrences).  Order p costs O(p^2) operations
 * per tape node.
 *
 * Used by the arbitrary-precision NDSolve driver (ndsolve_mpfr.c) as a
 * high-order integrator whose steps come with exact local derivatives -- the
 * data a high-order Hermite InterpolatingFunction needs. */
#ifndef MATHILDA_NDSOLVE_TAYLOR_H
#define MATHILDA_NDSOLVE_TAYLOR_H

#ifdef USE_MPFR
#include <stdbool.h>
#include <stddef.h>
#include <mpfr.h>
#include "../expr.h"

typedef struct NdTaylorTape NdTaylorTape;

/* Build the tape for the d right-hand sides f[0..d-1] over the state symbols
 * ysym[0..d-1] and the independent variable `tvar` (interned).  Constant
 * subexpressions are evaluated once at `bits`.  Returns NULL when some
 * subexpression is outside the supported vocabulary (the caller then uses a
 * Runge-Kutta method instead).  Inputs are borrowed. */
NdTaylorTape* nd_taylor_build(Expr* const* f, size_t d, Expr* const* ysym,
                              const char* tvar, long bits);

void nd_taylor_free(NdTaylorTape* T);

/* Taylor coefficients of the solution through (t0, Y0) up to order p:
 * C[k*d + i] receives Y_i^(k)(t0) / k! for k = 0..p (C must hold (p+1)*d
 * initialised mpfr_t).  Returns false when a coefficient is undefined (a pole
 * or a branch point of the RHS at the expansion point, or a non-finite
 * value). */
bool nd_taylor_coeffs(NdTaylorTape* T, const mpfr_t t0, mpfr_t* Y0, int p,
                      mpfr_t* C);

/* Machine-precision twin of nd_taylor_coeffs: C[k*d + i] = Y_i^(k)(t0)/k!
 * in double, k = 0..p.  Meant for the low orders a dense-output table needs
 * (y'', y''' at the nodes of a machine-precision solution). */
bool nd_taylor_coeffs_d(NdTaylorTape* T, double t0, const double* Y0, int p,
                        double* C);

#endif /* USE_MPFR */
#endif

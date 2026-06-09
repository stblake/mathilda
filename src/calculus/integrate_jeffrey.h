/* integrate_jeffrey.h
 *
 * Jeffrey-Rich continuous Weierstrass-substitution integrator for rational
 * functions of trigonometric (and, by generalisation, hyperbolic) functions.
 *
 * Reference: D. J. Jeffrey and A. D. Rich, "The Evaluation of Trigonometric
 * Integrals Avoiding Spurious Discontinuities", ACM Transactions on
 * Mathematical Software 20(1), 1994, pp. 124-135.
 *
 * The classical tan(x/2) substitution turns a rational function of sin x and
 * cos x into a rational function of u, which is integrated and back-substituted.
 * That step alone introduces spurious jump discontinuities at the poles of
 * tan(x/2) (odd multiples of pi).  This module removes them by adding the
 * Jeffrey-Rich secular correction  K * Floor[(x - b)/p], where K is the size of
 * the jump computed from the limits of the u-antiderivative at +/- infinity.
 *
 * Generalised to hyperbolic integrands via u = tanh(x/2); because tanh(x/2) is a
 * monotone bijection R -> (-1, 1) with no poles, that substitution introduces no
 * spurious discontinuity and needs no correction.
 *
 * Exposed as Integrate[f, x, Method -> "Weierstrass"] and Integrate`Weierstrass.
 */
#ifndef INTEGRATE_JEFFREY_H
#define INTEGRATE_JEFFREY_H

#include "expr.h"

/* Cascade entry (Automatic).  Fires only for genuine rational functions of the
 * trig/hyperbolic kernels of x that carry a kernel in a denominator, so the
 * cleaner table/Risch forms for polynomial trig (Integrate[Sin[x], x] etc.) are
 * left untouched.  Returns an owned antiderivative, or NULL when not applicable.
 * Borrows f and x. */
Expr* integrate_jeffrey_try(Expr* f, Expr* x);

/* Explicit-method entry (Method -> "Weierstrass"): same as the builtin but
 * called directly with (f, x).  No denominator gate.  Borrows f and x. */
Expr* integrate_jeffrey_full(Expr* f, Expr* x);

/* Explicit Method -> "Weierstrass" / Integrate`Weierstrass[f, x].  Applies to
 * any rational function of the trig/hyperbolic kernels of x (no denominator
 * gate).  Returns NULL when f is not of that form or the reduced rational
 * integral does not close. */
Expr* builtin_integrate_jeffrey(Expr* res);

/* Registration: installs Integrate`Weierstrass with attributes + docstring. */
void integrate_jeffrey_init(void);

#endif /* INTEGRATE_JEFFREY_H */

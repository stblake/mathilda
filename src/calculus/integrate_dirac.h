/*
 * integrate_dirac.h -- Sifting property of DiracDelta under a definite integral.
 *
 *   Integrate[DiracDelta[a x + b] h(x), {x, lo, hi}]
 *       = h(x0) / |a| * B(x0; lo, hi),          x0 = -b/a
 *
 * where B is the fraction of the impulse captured by the interval:
 *   - both limits numeric:  1 (lo < x0 < hi) / 0 (outside) / 1/2 (at an endpoint);
 *   - symbolic upper limit (the variation-of-parameters convolution case,
 *     lo numeric, hi = t, x0 >= lo):  HeavisideTheta[hi - x0], i.e. the FULL step
 *     even when x0 == lo -- the causal Green's-function convention that matches
 *     Mathematica's DSolve output for impulse forcing (e.g. x'' + 4 x == d(t)
 *     gives HeavisideTheta[t] Sin[2 t]/2).
 *
 * A sum integrand (impulse forcing mixed with an ordinary term, e.g. 1 + d(t-2)
 * or t + d(t)) is handled by ExpandAll + linearity: each DiracDelta term is
 * sifted, the delta-free remainder is re-integrated by the normal cascade.
 *
 * This is the only mechanism that assigns a value to a DiracDelta integrand;
 * DiracDelta is otherwise inert, so the recognizer can only affect currently
 * unevaluated inputs (no regression to any existing Integrate result).
 */

#ifndef MATHILDA_INTEGRATE_DIRAC_H
#define MATHILDA_INTEGRATE_DIRAC_H

#include "expr.h"
#include <stdbool.h>

/* Definite integral of an integrand carrying a DiracDelta[linear-in-x] factor.
 * Borrowed args (f the integrand, x the variable, a/b the lower/upper limits).
 * Returns a fresh value, or NULL to fall through -- including when f is
 * delta-free, when the delta argument is nonlinear, or when the interval
 * membership of the impulse cannot be decided. */
Expr* integrate_dirac_try(Expr* f, Expr* x, Expr* a, Expr* b);

/* Does f carry a DiracDelta whose argument depends on x? When it does and
 * integrate_dirac_try declines (a multiple root such as DiracDelta[x^2], two
 * impulses in x, an undecided endpoint), the definite integral stays
 * unevaluated: the other methods read the inert delta as an ordinary function
 * that is 0 almost everywhere and answered 0. */
bool integrate_dirac_in(Expr* f, Expr* x);

#endif /* MATHILDA_INTEGRATE_DIRAC_H */

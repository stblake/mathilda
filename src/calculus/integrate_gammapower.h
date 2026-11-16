/* integrate_gammapower.h — symbolic power times an exponential -> incomplete Gamma.
 *
 * K x^p E^(a x^m) with the exponent p SYMBOLIC (not a number) -> the incomplete
 * Gamma[s, z]; and (v0.6588) K x^p T(c x^q), T in {Exp, Sin, Cos, Sinh, Cosh},
 * q an exact nonzero rational, s = (p+1)/q with 2s not an integer, through the
 * branch-exact -(1/q) x^(p+1) w^(-s) Gamma[s, w], w = -a x^q.  Returns a fresh,
 * diff-back-verified antiderivative or NULL.
 * Defined in integrate_gammapower.c.
 */

#ifndef MATHILDA_INTEGRATE_GAMMAPOWER_H
#define MATHILDA_INTEGRATE_GAMMAPOWER_H

#include "expr.h"

Expr* integrate_gammapower_try(Expr* f, Expr* x);

#endif /* MATHILDA_INTEGRATE_GAMMAPOWER_H */

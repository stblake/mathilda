/* integrate_gammapower.h — symbolic power times an exponential -> incomplete Gamma.
 *
 * K x^p E^(a x^m) with the exponent p SYMBOLIC (not a number) -> the incomplete
 * Gamma[s, z].  Returns a fresh, diff-back-verified antiderivative or NULL.
 * Defined in integrate_gammapower.c.
 */

#ifndef MATHILDA_INTEGRATE_GAMMAPOWER_H
#define MATHILDA_INTEGRATE_GAMMAPOWER_H

#include "expr.h"

Expr* integrate_gammapower_try(Expr* f, Expr* x);

#endif /* MATHILDA_INTEGRATE_GAMMAPOWER_H */

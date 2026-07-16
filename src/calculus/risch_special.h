/* risch_special.h — special-function outputs & base cases for the Risch integrator.
 *
 * The non-tower cases that emit higher transcendental functions or delegate to
 * the rational base: the pure-rational case (delegated to
 * Integrate`BronsteinRational), the Erf/Erfi, ExpIntegralEi, LogIntegral, and
 * dilogarithm (PolyLog) recognisers, and the logarithmic-polynomial case.
 * These are tried before the tower/exponential machinery.  Defined in
 * risch_special.c.
 */

#ifndef MATHILDA_RISCH_SPECIAL_H
#define MATHILDA_RISCH_SPECIAL_H

#include "expr.h"

/* Pure rational function of x -> Bronstein rational integrator (NULL if not
 * rational in x). */
Expr* rt_rational_case(Expr* f, Expr* x);

/* Erf / Erfi / ExpIntegralEi / LogIntegral / dilogarithm recogniser. */
Expr* rt_special_case(Expr* f, Expr* x);

/* Logarithmic-polynomial case (polynomial in a single Log[u]). */
Expr* rt_log_poly_case(Expr* f, Expr* x);

#endif /* MATHILDA_RISCH_SPECIAL_H */

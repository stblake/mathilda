/* integrate_logbyparts.h — Log[g] times a trig/hyperbolic kernel -> one
 * integration by parts (Si/Ci/… residual).
 *
 * Recognizes  c Log[g(x)] K(x)  where K carries a trig or hyperbolic kernel of x
 * and is free of Log, and integrates it by one step of integration by parts,
 *
 *     INT Log[g] K dx  =  Log[g] V  -  INT V (g'/g) dx,    V = INT K dx,
 *
 * so the Log[x] Sin[a x] / Log[x] Cos[a x] family closes to the clean
 * Log*trig + SinIntegral/CosIntegral form instead of falling through to the
 * seconds-long ParallelMixedSpecial search.  Returns a fresh, diff-back-verified
 * antiderivative or NULL (decline).  Defined in integrate_logbyparts.c.
 */

#ifndef MATHILDA_INTEGRATE_LOGBYPARTS_H
#define MATHILDA_INTEGRATE_LOGBYPARTS_H

#include "expr.h"

Expr* integrate_logbyparts_try(Expr* f, Expr* x);

/* c S K with exactly one S in {ExpIntegralEi, Erf, Erfi, Erfc, LogIntegral} of
 * an x-dependent argument and K free of special functions: INT S K = S V -
 * INT V S', V = INT K elementary (INT E^x Ei[x] = E^x Ei[x] - Ei[2 x]).
 * PossibleZeroQ-verified, or NULL. */
Expr* integrate_specialbyparts_try(Expr* f, Expr* x);

#endif /* MATHILDA_INTEGRATE_LOGBYPARTS_H */

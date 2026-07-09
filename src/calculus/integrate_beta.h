/*
 * integrate_beta.h -- Definite integration by Euler-Beta reduction.
 *
 * Two related closed-form families that the elementary (Newton-Leibniz)
 * cascade cannot reach because the antiderivatives are non-elementary (an
 * incomplete Beta), yet the *definite* value over the canonical interval is a
 * Beta function:
 *
 *   Beta on [0,1]   Integrate[x^(k-1) (1-x)^(l-1), {x,0,1}] = Beta[k, l],
 *                   with Log[x]^i Log[1-x]^j weights -> the mixed parameter
 *                   derivative d^i/dk^i d^j/dl^j Beta[k, l].
 *
 *   Trig powers     Integrate[Sin[x]^m Cos[x]^n, {x,0,Pi/2}]
 *                     = Beta[(m+1)/2, (n+1)/2] / 2,
 *                   with the [0,Pi] and [0,2Pi] cases handled by the standard
 *                   parity multipliers (odd power -> 0, else 2x / 4x the
 *                   quarter-period value).
 *
 * Both are gated on the Beta convergence strip (Re of each argument > 0); when
 * Assumptions do not settle it the result is a ConditionalExpression.  Integer
 * powers already close through Newton-Leibniz (their antiderivative IS
 * elementary), so these mechanisms run after it and mainly serve the
 * non-integer / symbolic / even-over-[0,Pi] cases.
 *
 * Reachable via Method -> "Beta" / "TrigPower" and the explicit entry points
 *   Integrate`Beta[f, {x, 0, 1}]
 *   Integrate`TrigPower[f, {x, 0, c}]
 */

#ifndef MATHILDA_INTEGRATE_BETA_H
#define MATHILDA_INTEGRATE_BETA_H

#include "expr.h"

/* Beta[k,l] family on [0,1] (with Log[x]/Log[1-x] weights).  Borrowed args;
 * assumptions may be NULL.  Returns a fresh value or NULL to fall through. */
Expr* integrate_beta_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);

/* Sin^m Cos^n family over [0, Pi/2], [0, Pi], [0, 2 Pi].  Borrowed args;
 * assumptions may be NULL.  Returns a fresh value or NULL to fall through. */
Expr* integrate_trigpower_try(Expr* f, Expr* x, Expr* a, Expr* b,
                              Expr* assumptions);

/* Explicit builtins.  Strict: NULL on any non-applicable input. */
Expr* builtin_integrate_beta(Expr* res);
Expr* builtin_integrate_trigpower(Expr* res);

/* Register both builtins + attributes + docstrings. */
void integrate_beta_init(void);

#endif /* MATHILDA_INTEGRATE_BETA_H */

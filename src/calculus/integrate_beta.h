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

/* Euler integral x^(a-1)(1-x)^(b-1)(alpha+beta x)^e on [0,1] -> the Euler
 * representation alpha^e Beta[a,b] Hypergeometric2F1[-e, a, a+b, -beta/alpha].
 * Requires the extra linear factor; pure Beta is integrate_beta_try's.  Borrowed
 * args; assumptions may be NULL.  Returns a fresh value or NULL to fall through. */
Expr* integrate_euler_2f1_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);

/* Sin^m Cos^n family over [0, Pi/2], [0, Pi], [0, 2 Pi].  Borrowed args;
 * assumptions may be NULL.  Returns a fresh value or NULL to fall through. */
Expr* integrate_trigpower_try(Expr* f, Expr* x, Expr* a, Expr* b,
                              Expr* assumptions);

/* Log[Sin[x]] / Log[Cos[x]] over [0, Pi/2] (both -> -(Pi/2) Log 2) or
 * [0, Pi] (Log[Sin] -> -Pi Log 2; Log[Cos] complex, declined).  The trig
 * argument must be exactly x.  Borrowed args; returns fresh or NULL. */
Expr* integrate_logtrig_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);

/* Poisson log integral Integrate[Log[A + B Cos[x]], {x,0,Pi}] =
 * Pi Log[(A + Sqrt[A^2-B^2])/2] (A>0, A>=|B|; A,B free of x).  Borrowed args. */
Expr* integrate_poissonlog_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);

/* King's-property reflection for Integrate[x g(x), {x,0,c}] with g invariant
 * under x -> c-x:  = (c/2) Integrate[g, {x,0,c}].  Borrowed args. */
Expr* integrate_kingreflect_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);

/* Raabe's integral: Integrate[Log[Gamma[x]], {x, t, t+1}] = (1/2)Log[2 Pi]
 * + t Log[t] - t for t >= 0 (so [0,1] gives (1/2)Log[2 Pi]). */
Expr* integrate_raabe_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);

/* Log-weighted rational on [0,1]: Integrate[x^(c-1) Log[x]^n/(1-x^m), {x,0,1}]
 * = (-1)^n n!/m^(n+1) HurwitzZeta[n+1, c/m]  (c>0, m>=1, n>=1). */
Expr* integrate_ratlog01_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);

/* Half-Beta on [0,1]: Integrate[C x^(c-1)(1-x^2)^e Log[x]^i, {x,0,1}]
 * = C (1/2)^(i+1) d^i/dA^i Beta[A, e+1]|_{A=c/2} (via t = x^2). */
Expr* integrate_halfbeta_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);

/* t = Tan[x] on [0, b], 0 < b <= Pi/4, for a trig integrand with a symbolic
 * exponent, handed on to the [0,1] Beta stages above.  Borrowed args; returns
 * a fresh value or NULL to fall through. */
Expr* integrate_tanbeta_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);
/* By parts on [0,1] for C x^p ArcSin[x] / ArcCos[x], p symbolic: boundary term
 * plus a Half-Beta remainder.  Borrowed args; fresh value or NULL. */
Expr* integrate_invtrig_beta_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);
/* Explicit builtins.  Strict: NULL on any non-applicable input. */
Expr* builtin_integrate_beta(Expr* res);
Expr* builtin_integrate_trigpower(Expr* res);

/* Register both builtins + attributes + docstrings. */
void integrate_beta_init(void);

#endif /* MATHILDA_INTEGRATE_BETA_H */

/*
 * integrate_elliptic.h -- definite integrals of R(x, Sqrt[P(x)]) with P a
 * squarefree cubic or quartic, in closed form through the Legendre elliptic
 * integrals EllipticF / EllipticE / EllipticPi / EllipticK plus elementary
 * terms.
 *
 * Method (see integrate_elliptic.c for the derivations):
 *   1. normalise the integrand to A(x) + W(x)/Sqrt[P(x)] on the interval
 *      (every half-integer power of a rational function sharing one
 *      squarefree kernel P is folded onto a single Sqrt[P], P > 0 there);
 *   2. a real Moebius map x = (al t + be)/(ga t + de) sending 0 and Infinity
 *      to the common harmonic pair of the two real quadratic factors of P
 *      makes P even in t: (ga t + de)^4 P = pi (A1 + B1 t^2)(A2 + B2 t^2);
 *      this covers four real roots, two real roots and a complex pair, two
 *      complex pairs and (with Infinity as the fourth root) every cubic;
 *   3. Hermite reduction in t leaves the basis Int dt/Sqrt, Int t dt/Sqrt,
 *      Int t^2 dt/Sqrt and Int dt/((t - tau) Sqrt) plus an algebraic part;
 *   4. s = t^2 turns the even parts into a cubic with three REAL linear
 *      factors s, A1 + B1 s, A2 + B2 s and the odd parts into a quadratic
 *      radical, evaluated with Carlson's reduction formulas (DLMF 19.29.4,
 *      19.29.7, 19.29.8; Cauchy principal values where the even/odd split
 *      creates a removable pole inside the range), and converted to
 *      Legendre form exactly (DLMF 19.25.5-19.25.14, 19.20.14, 19.2.17-20).
 *   5. the closed form is checked at 30 digits against an independent
 *      90-digit double-exponential quadrature of the integrand; any
 *      disagreement declines (never a wrong answer).
 *
 * Reachable as Integrate`Elliptic[f, {x, a, b}] (optionally with
 * Assumptions -> ...); returns unevaluated when it cannot evaluate.
 */

#ifndef MATHILDA_INTEGRATE_ELLIPTIC_H
#define MATHILDA_INTEGRATE_ELLIPTIC_H

#include "expr.h"

/* Core entry point.  f, x, a, b, assumptions are borrowed; assumptions may
 * be NULL.  Returns a fresh closed form, or NULL to decline. */
Expr* integrate_elliptic_try(Expr* f, Expr* x, Expr* a, Expr* b,
                             Expr* assumptions);

/* Integrate`Elliptic[f, {x, a, b}, opts...] builtin. */
Expr* builtin_integrate_elliptic(Expr* res);

/* Register the builtin, its attributes and docstring. */
void integrate_elliptic_init(void);

/* Int over [0, Pi/2] (or [0, Pi], [0, 2 Pi]) of (A + B Sin[x]^2)^(+-1/2), A > 0,
 * as EllipticK / EllipticE; NULL when the integrand is not of that shape. */
Expr* integrate_elliptic_trig_quadratic(Expr* f, Expr* x, Expr* a, Expr* b);

#endif /* MATHILDA_INTEGRATE_ELLIPTIC_H */

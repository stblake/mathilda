/*
 * integrate_intrep.h -- Definite integration by recognising a classical
 * INTEGRAL REPRESENTATION of a special function.
 *
 * A family of half-line integrals is, by a standard identity, the DEFINITION of
 * a special function on its convergence domain; the integrand has no elementary
 * antiderivative, so neither Newton-Leibniz nor the Mellin/residue families
 * close it, but the closed form is immediate once the shape is matched and the
 * convergence gate is proved from the Assumptions.  Correct-by-construction (the
 * representation is an identity on its domain); NO NIntegrate crosscheck.
 *
 *   Laplace-Bessel   Integrate[E^(-p x) BesselJ[nu, q x], {x,0,Inf}]
 *                      = (Sqrt[q^2+p^2]-p)^nu / (q^nu Sqrt[q^2+p^2]),   p>0.
 *   Bessel-K (cosh)  Integrate[E^(-A Cosh[x]) Cosh[n x], {x,0,Inf}]
 *                      = BesselK[n, A],                                 A>0.
 *   Bessel-K (exp)   Integrate[x^(nu-1) E^(-A x - B/x), {x,0,Inf}]
 *                      = 2 (B/A)^(nu/2) BesselK[nu, 2 Sqrt[A B]],       A>0,B>0.
 *   Airy             Integrate[Cos[p x^3 + q x], {x,0,Inf}]
 *                      = Pi (3p)^(-1/3) AiryAi[q (3p)^(-1/3)],          p>0, q real.
 *
 * Reachable via Method -> "IntegralRepresentation" and the explicit entry point
 *   Integrate`IntegralRepresentation[f, {x, 0, Infinity}]
 */

#ifndef MATHILDA_INTEGRATE_INTREP_H
#define MATHILDA_INTEGRATE_INTREP_H

#include "expr.h"

/* Integral-representation recognizer on the half line [0, Infinity).  Borrowed
 * args; assumptions may be NULL.  Returns a fresh value or NULL to fall through. */
Expr* integrate_intrep_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);

/* Explicit builtin.  Strict: NULL on any non-applicable input. */
Expr* builtin_integrate_intrep(Expr* res);

/* Register the builtin + attributes + docstring. */
void integrate_intrep_init(void);

#endif /* MATHILDA_INTEGRATE_INTREP_H */

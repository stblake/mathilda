/* integrate_fresnel.h — Fresnel integrals of a Gaussian-phase trig integrand.
 *
 * K Sin[a x^2 + b x + c] and K Cos[...] (a != 0, K free of x) -> FresnelS/FresnelC
 * by completing the square, the trigonometric sibling of the K E^(a x^2+...) ->
 * Erf recognizer.  Returns a fresh, diff-back-verified antiderivative or NULL.
 * Defined in integrate_fresnel.c.
 */

#ifndef MATHILDA_INTEGRATE_FRESNEL_H
#define MATHILDA_INTEGRATE_FRESNEL_H

#include "expr.h"

Expr* integrate_fresnel_try(Expr* f, Expr* x);

#endif /* MATHILDA_INTEGRATE_FRESNEL_H */

/* integrate_linarg.h
 *
 * Linear-argument substitution stage for the indefinite-integration cascade.
 *
 * Recognises an integrand that is a rational function of trig / hyperbolic
 * kernels all sharing ONE non-trivial linear argument w = a*x + b (a a non-zero
 * number, and a != 1 or b != 0), with a kernel in a denominator, and reduces it
 * to the bare-argument integral via u = a*x + b:
 *
 *     Integrate[f(a*x+b), x]  ->  (1/a) * ( Integrate[f(u), u] /. u -> a*x+b ).
 *
 * The substitution is an exact identity (dx = du/a is a constant Jacobian), so
 * the result is correct by construction whenever the bare sub-integral closes.
 *
 * Why it exists: the Weierstrass stage (integrate_jeffrey.c) always substitutes
 * t = Tan[x/2], so a scaled argument such as Sec[3x]^2 is first multiple-angle
 * expanded into a high-degree rational in Tan[x/2] and comes back a degree-12
 * mess (and Integrate[Sec[3x]^2, x] = Tan[3x]/3 would then never appear). This
 * stage runs IMMEDIATELY BEFORE Weierstrass: it rewrites the scaled/shifted
 * argument to a bare one, which the recursive Integrate then closes cleanly
 * (reusing the identical bare-argument machinery), and declines otherwise so
 * Weierstrass still gets its turn. Polynomial trig (Sin[2x]^3) and bare-argument
 * integrands are left untouched (no kernel in a denominator / trivial w).
 *
 * Returns a new owned Expr* on success, or NULL to continue the cascade.
 */
#ifndef INTEGRATE_LINARG_H
#define INTEGRATE_LINARG_H

#include "expr.h"

Expr* integrate_linarg_try(Expr* f, Expr* x);

#endif /* INTEGRATE_LINARG_H */

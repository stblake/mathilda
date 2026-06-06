/*
 * integrate_interp.h
 *
 * Indefinite integration of an applied InterpolatingFunction object:
 *
 *     Integrate[InterpolatingFunction[...][x], x]
 *
 * mirrors the way D[InterpolatingFunction[...][x], x] differentiates such
 * objects.  Whereas differentiation only bumps the derivative-order
 * annotation, integration must build a genuinely new interpolant — the
 * antiderivative — because the per-window evaluation kernels cannot evaluate
 * negative ("anti-") derivative orders.
 */
#ifndef INTEGRATE_INTERP_H
#define INTEGRATE_INTERP_H

#include "../expr.h"

/*
 * integrate_interp:
 *   `f` is the integrand and `x` the (symbol) integration variable, exactly as
 *   handed to builtin_integrate.  If `f` is an applied 1-D InterpolatingFunction
 *   object whose argument is `x` itself, returns the antiderivative as a fresh
 *   `InterpolatingFunction[...][x]` expression.  Returns NULL for anything it
 *   does not handle (non-ifun integrand, multidimensional object, an applied
 *   argument other than the bare variable, …), leaving the caller to continue
 *   its normal integration cascade.  Both inputs are borrowed.
 */
Expr* integrate_interp(Expr* f, Expr* x);

#endif /* INTEGRATE_INTERP_H */

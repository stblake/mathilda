/*
 * interp.h
 *
 * InterpolatingFunction: an approximate-function object whose values are
 * found by piecewise-polynomial interpolation of tabulated data.
 *
 *   InterpolatingFunction[domain, table]
 *
 *     domain = {{xmin, xmax}}            -- the interval the data span.
 *     table  = {{x1, y1}, ..., {xn, yn}} -- data points, xi strictly increasing.
 *
 * The object itself is a normal form (like Function[...]) and persists
 * unevaluated. Applying it, InterpolatingFunction[...][x], evaluates the
 * approximate function at x via interp_apply().
 */
#ifndef INTERP_H
#define INTERP_H

#include "expr.h"

/* Register InterpolatingFunction in the symbol table (attributes only). */
void interp_init(void);

/*
 * interp_apply:
 *   `ifun` is the InterpolatingFunction[...] object (an EXPR_FUNCTION whose
 *   head symbol is InterpolatingFunction). `call_args`/`argc` are the
 *   arguments it is being applied to, e.g. the `x` in ifun[x].
 *
 *   Returns a freshly-allocated Expr* holding the interpolated value, or
 *   NULL if the application cannot be evaluated (non-numeric argument,
 *   malformed object, wrong arity) — in which case the caller leaves the
 *   expression unevaluated. The caller retains ownership of every input.
 */
Expr* interp_apply(Expr* ifun, Expr** call_args, size_t argc);

/*
 * interp_make_derivative:
 *   Reduce Derivative[d1, ..., dm][InterpolatingFunction[...]] to a fresh
 *   InterpolatingFunction object carrying the accumulated derivative orders
 *   (additive with any orders the object already holds). `deriv_head` is the
 *   Derivative[d1, ..., dm] head; `ifun` is the InterpolatingFunction object.
 *   Returns a newly-allocated object, or NULL when the order count does not
 *   match the object's dimensionality. Both inputs are borrowed.
 */
Expr* interp_make_derivative(Expr* deriv_head, Expr* ifun);

#endif /* INTERP_H */

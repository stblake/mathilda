#ifndef SERIES_H
#define SERIES_H

#include "expr.h"

/*
 * Registers SeriesData, Series, and Normal.
 *
 *   SeriesData[x, x0, {a0, a1, ...}, nmin, nmax, den]
 *     Represents a power series in x about the point x0. The coefficient
 *     ai is attached to (x - x0)^((nmin + i)/den), and an
 *     O[x - x0]^(nmax/den) term represents the dropped higher-order terms.
 *     SeriesData is a pure data head with attribute Protected.
 *
 *   Series[f, {x, x0, n}]
 *     Produces a truncated power-series expansion of f about x = x0 up to
 *     order (x - x0)^n. Handles Taylor, Laurent (negative powers), Puiseux
 *     (fractional powers), and logarithmic expansions. Attributes:
 *     HoldAll, Protected.
 *
 *   Series[f, x -> x0]
 *     Emits the leading term only (expansion at order 0).
 *
 *   Series[f, {x, x0, nx}, {y, y0, ny}, ...]
 *     Iterated multivariate expansion, right-to-left.
 *
 *   Normal[SeriesData[...]]
 *     Drops the O-term and returns the ordinary sum expression.
 */
void series_init(void);

Expr* builtin_series(Expr* res);
Expr* builtin_normal(Expr* res);

/*
 * series_split_two_term
 *   Decompose e = a + b*x^(exp_num/exp_den) structurally, without running
 *   the full series-expansion pipeline. Returns true on success and
 *   transfers ownership of *a_out / *b_out to the caller. Returns false
 *   and leaves the out-pointers NULL otherwise (the shape didn't match).
 *
 * Handles:
 *   - e == x                       -> (0, 1, 1/1)
 *   - expressions free of x        -> (e, 0, 1/1)
 *   - Plus[summands]               -> sum of a-parts; b-parts must share one exponent
 *   - Times[factors]               -> at most one non-free factor
 *   - Power[x, rational]           -> (0, 1, p/q)
 *
 * Exposed for unit testing of the probe itself; the caller never owns `x`.
 */
bool series_split_two_term(Expr* e, Expr* x,
                           Expr** a_out, Expr** b_out,
                           int64_t* exp_num, int64_t* exp_den);

/*
 * series_differentiate / series_integrate
 *   Differentiate or integrate a SeriesData[x, x0, {coefs}, nmin, nmax, den]
 *   object term-by-term in `var`. When `var` is the series variable the power
 *   rule is applied to each term (integration uses constant 0, as in
 *   Mathematica); when `var` is a different symbol the coefficients are
 *   transformed and the powers kept (valid only when x0 is free of `var`).
 *   Return a fresh SeriesData Expr on success, or NULL when the shape is
 *   unsupported (e.g. non-integer orders) or, for integration, when a nonzero
 *   (x-x0)^-1 term would require a Log. The caller then leaves D[]/Integrate[]
 *   unevaluated.
 */
Expr* series_differentiate(Expr* sd, Expr* var);
Expr* series_integrate(Expr* sd, Expr* var);

/*
 * SeriesData arithmetic
 *   These drive the internal SeriesObj algebra (so_add / so_mul / so_pow_*)
 *   from the arithmetic builtins so that Plus, Times and Power combine power
 *   series the way Mathematica does. Divide and Subtract are covered for free
 *   because they rewrite to Times[a, Power[b,-1]] and Plus[a, Times[-1,b]].
 *
 * is_series_data
 *   Cheap structural test: e is SeriesData[x, x0, {coefs}, nmin, nmax, den].
 *
 * series_combine_plus / series_combine_times
 *   Fold the operands of a Plus / Times node (at least one of which is a
 *   SeriesData) into a single SeriesData. Operands that are not series are
 *   series-expanded about the common (x, x0) to the controlling (minimum)
 *   order. Returns a fresh SeriesData Expr, or NULL when the operands are
 *   incompatible (different expansion variable or point, or a non-series
 *   operand that cannot be expanded) -- the caller then leaves the node
 *   unevaluated, matching Mathematica's behaviour for e.g. two series about
 *   different points.
 *
 * series_power
 *   base^exp where base and/or exp is a SeriesData. Integer exponents (incl.
 *   negative, via series inversion) and exponents free of the series variable
 *   are handled directly; an exponent that depends on the series variable (or
 *   a non-series base raised to a series exponent, e.g. 2^series) is rewritten
 *   as Exp[exp*Log[base]] and re-expanded. Returns NULL when the shape is
 *   unsupported.
 */
bool  is_series_data(const Expr* e);
Expr* series_combine_plus(Expr* const* args, size_t n);
Expr* series_combine_times(Expr* const* args, size_t n);
Expr* series_power(Expr* base, Expr* exp);

#endif // SERIES_H

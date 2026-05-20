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

#endif // SERIES_H

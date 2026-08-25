/*
 * reduce_real_util.h
 *
 * Shared real-algebraic primitives for `Reduce`'s Reals engines (the univariate
 * sign diagram, reduce_univar.c, and the multivariate CAD, reduce_cad.c).
 *
 * Everything here is exact: numeric literals are signed directly and every other
 * constant real value goes through `flint_qqbar_compare` (the exact real-
 * algebraic oracle), which returns -2 when it cannot decide (non-real, non-
 * algebraic, or FLINT unavailable).  A -2 must propagate to a NULL / decline in
 * the caller -- the engines never guess a sign.
 *
 * Root isolation is delegated to `Solve[..., Reals]` (rationals, radicals, and
 * `Root[...]` objects), and rational sample points are certified with the sign
 * oracle so a rational-root problem needs no FLINT.
 */
#ifndef REDUCE_REAL_UTIL_H
#define REDUCE_REAL_UTIL_H

#include "expr.h"
#include "reduce_opts.h"
#include <stdbool.h>

/* Exact sign of a constant real value: -1, 0, 1, or -2 (undecidable / FLINT
 * off).  Plain numbers (Integer/Real/BigInt/MPFR/Rational) are signed directly;
 * anything else goes to the qqbar oracle. */
int rru_sign_of(const Expr* e);

/* Sign of (a - b): -1, 0, 1, or -2.  Computed by evaluating a - b then signing
 * the constant result. */
int rru_sign_compare(const Expr* a, const Expr* b);

/* Numeric double approximation of `e` via N[]; sets *ok=false when the result is
 * not a machine-representable number. */
double rru_approx_double(const Expr* e, bool* ok);

/* A rational strictly between `lo` and `hi` (either NULL means unbounded on that
 * side), certified exactly with rru_sign_compare.  Returns an owned Integer or
 * Rational Expr, or NULL if none is found within the search budget or a
 * comparison is undecidable. */
Expr* rru_rational_between(const Expr* lo, const Expr* hi);

/* True iff `poly` is a polynomial in `x` (PolynomialQ[poly, x]). */
bool rru_is_polynomial(const Expr* poly, const Expr* x);

/* Sign of poly with the single rule x -> sample applied: -1, 0, 1, or -2. */
int rru_poly_sign_at(const Expr* poly, const Expr* x, const Expr* sample);

/* Emit the union of satisfying cells of a 1-D sign diagram as a logical Expr.
 * `bp` holds the m ascending breakpoints; `truth` has 2m+1 entries, one per cell
 * in geometric order (interval I0, point P0, interval I1, ..., point P_{m-1},
 * interval Im), each 1 (in the solution set) or 0.  Produces True / False / an
 * `x REL b` / a two-sided `Inequality` chain / an `Or` of segments / a cofinite
 * `And` of `x != b`, then one evaluate() pass to flatten.  Does not free `bp` or
 * `truth` (the caller owns them).  Shared by the polynomial sign diagram
 * (reduce_univar.c) and the elementary-real-function diagram (reduce_realdiag.c).
 */
Expr* rru_emit_sign_diagram(Expr** bp, int m, const int* truth, const Expr* x);

/* Append the real roots of `poly` (isolated via Solve over Reals) to *arr,
 * growing (*arr,*n,*cap).  Returns false to bail (an unexpected Solve shape or a
 * ConditionalExpression / parametric root).
 *
 * Provenance: when `prov` is non-NULL it is grown in lockstep with *arr and each
 * pushed root records `factor_id` -- so a caller isolating roots factor-by-factor
 * can later recover which factor produced each breakpoint (needed by CAD's
 * symbolic sector emission).  Pass prov=NULL, factor_id=0 when provenance is not
 * wanted.
 *
 * `opts` forwards the Cubics / Quartics radical flags onto the internal
 * `Solve[..., Reals]` root isolation (pass NULL to keep Solve's defaults). */
bool rru_collect_roots(const Expr* poly, const Expr* x,
                       Expr*** arr, int* n, int* cap,
                       int** prov, int factor_id,
                       const ReduceOpts* opts);

#endif /* REDUCE_REAL_UTIL_H */

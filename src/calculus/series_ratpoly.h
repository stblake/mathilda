/*
 * series_ratpoly.h -- truncated power-series arithmetic over Q for the series
 * engine (series.c), on FLINT's fmpq_poly.
 *
 * The series engine keeps a coefficient array of Expr* per series and does
 * every coefficient operation through the evaluator.  When every coefficient
 * is a rational number -- the series of a polynomial, of a rational function
 * with rational coefficients, of a composition of the elementary functions at
 * 0 -- that is O(n^2) evaluator calls per product where one truncated
 * polynomial multiply (mullow) does the same work in quasi-linear time.  These
 * entry points take and return plain coefficient arrays (index i = the i-th
 * coefficient from the leading exponent) so series.c needs no FLINT types.
 *
 * Every function returns 1 on success, filling `out[0..n-1]` with fresh owned
 * Exprs (Integer / BigInt / Rational), and 0 -- touching nothing -- when some
 * input coefficient is not a rational number, so the caller falls through to
 * its general path.  Without USE_FLINT they all return 0.
 */
#ifndef SERIES_RATPOLY_H
#define SERIES_RATPOLY_H

#include <stddef.h>
#include <stdint.h>
#include "expr.h"

/* True iff e is an Integer, BigInt or Rational. */
int srp_is_rational(const Expr* e);

/* out = (a * b) mod t^n, for a of length na and b of length nb. */
int srp_mullow(Expr* const* a, size_t na, Expr* const* b, size_t nb,
               size_t n, Expr** out);

/* out = 1/a mod t^n.  Requires a[0] != 0 (returns 0 otherwise). */
int srp_inv_series(Expr* const* a, size_t na, size_t n, Expr** out);

/* out = a^e mod t^n, e >= 0. */
int srp_pow_trunc(Expr* const* a, size_t na, uint64_t e, size_t n, Expr** out);

/* out = K(u) mod t^n where K = sum_k k[k] s^k (nk terms) and u has u[0] == 0
 * (the series u as coefficients of t^0 .. t^(nu-1)); returns 0 if u[0] != 0. */
int srp_compose_series(Expr* const* k, size_t nk, Expr* const* u, size_t nu,
                       size_t n, Expr** out);

/* a + b and a * b of two rational numbers, as fresh owned Exprs; NULL when
 * either is not a rational number. */
Expr* srp_add2(const Expr* a, const Expr* b);
Expr* srp_mul2(const Expr* a, const Expr* b);

#endif /* SERIES_RATPOLY_H */

#ifndef INTERVAL_H
#define INTERVAL_H

/* Interval arithmetic — Interval[{min, max}] and the union form
 * Interval[{a1,b1}, {a2,b2}, ...].
 *
 * An Interval is an ordinary EXPR_FUNCTION with head SYM_Interval whose
 * arguments are two-element List pairs, recognized by is_interval() and built
 * by make_interval_pair() — the same "numeric container is a function head"
 * pattern used for Complex[re,im] and Rational[n,d] (see src/arithmetic.c).
 *
 * Exact endpoints (Integer/BigInt/Rational, symbolic values such as Sin[2],
 * and +/-Infinity) are kept exact/symbolic; inexact endpoints (Real/MPFR) are
 * rounded OUTWARD (lower bound down, upper bound up) so a computed interval is
 * always a rigorous enclosure. Outward rounding is done by evaluating the
 * endpoint operation with the ordinary (correctly-rounded) kernels and then
 * nudging the inexact result one ULP outward — see iv_widen() in interval.c.
 */

#include "expr.h"
#include <stdbool.h>
#include <stdint.h>

/* Monotone direction for interval_thread_monotone. */
typedef enum { IV_INC, IV_DEC } IvMonotoneDir;

/* True iff e is Interval[...] with at least one argument. */
bool is_interval(const Expr* e);

/* Build Interval[{lo, hi}] taking ownership of lo and hi. */
Expr* make_interval_pair(Expr* lo, Expr* hi);

/* Number of disjoint pairs in a (canonical) interval. */
size_t interval_pair_count(const Expr* iv);
/* Borrowed lower / upper endpoint of the k-th pair (0-based). */
Expr* interval_pair_lo(const Expr* iv, size_t k);
Expr* interval_pair_hi(const Expr* iv, size_t k);

/* Registered builtin: canonicalizes Interval[...] (orders endpoints, sorts and
 * merges overlapping/touching pairs). Returns NULL when nothing changes. */
Expr* builtin_interval(Expr* res);

/* Module registration (symtab_add_builtin + attributes + docstring). */
void interval_init(void);

/* Compare two endpoint expressions numerically.
 * Returns the sign of (a - b): -1, 0, or +1, and sets *decidable to whether the
 * comparison could be resolved. Handles +/-Infinity, exact numbers, and
 * symbolic-numeric endpoints (Pi, Sin[2], ...) via a numeric fallback. */
int interval_endpoint_cmp(Expr* a, Expr* b, bool* decidable);

/* Interval kernels. Each takes borrowed operands and returns a fresh, canonical
 * Interval (owned by the caller), or NULL to leave the operation symbolic. */
Expr* interval_add(const Expr* A, const Expr* B);
Expr* interval_negate(const Expr* A);
Expr* interval_multiply(const Expr* A, const Expr* B);
Expr* interval_reciprocal(const Expr* A);
Expr* interval_power_int(const Expr* A, int64_t n);

/* A^p for a positive real/rational scalar exponent p when A is entirely >= 0
 * (e.g. Sqrt, which the parser lowers to Power[., 1/2]). Returns NULL otherwise. */
Expr* interval_power_pos_exp(const Expr* A, const Expr* p);

/* Thread a monotone unary function (given by head name, e.g. "Exp") over an
 * interval. Returns NULL if a result endpoint leaves the reals (e.g. Log of a
 * negative endpoint). */
Expr* interval_thread_monotone(const Expr* iv, IvMonotoneDir dir, const char* head);

/* Thread an elementary function over an interval, dispatching by head name:
 * monotone heads (Exp, Log, Sqrt, ArcTan, Sinh, Tanh) and the non-monotone
 * Sin/Cos/Tan/Abs (with critical-point range analysis). Returns NULL for heads
 * not yet supported, so callers fall through to symbolic. */
Expr* interval_apply_function(const char* head, const Expr* iv);

/* Thread PolyGamma[n, .] over an interval (monotone on (0, inf): increasing for
 * even n, decreasing for odd n). Returns NULL outside (0, inf) or for n < 0. */
Expr* interval_polygamma(int64_t n, const Expr* iv);

/* Promote a scalar number to the point interval Interval[{x, x}] (copies x). */
Expr* interval_from_scalar(const Expr* x);

/* Compare an interval against a scalar c for the comparison builtins.
 * Returns the sign of (interval - c): -1 if the interval lies entirely below c,
 * +1 if entirely above, and sets *decidable=false (returning 0) when c lies in
 * or touches the interval. */
int interval_compare_scalar(const Expr* iv, Expr* c, bool* decidable);

/* Compare two intervals. Returns the sign of (A - B): -1 if A lies entirely
 * below B, +1 if entirely above; 0 with *decidable=true only when A and B are
 * the same single point; *decidable=false (returning 0) otherwise. */
int interval_compare_intervals(const Expr* A, const Expr* B, bool* decidable);

/* base^Interval for a positive scalar base (monotone in the exponent). */
Expr* interval_scalar_pow_interval(const Expr* base, const Expr* B);
/* Interval^Interval for a strictly-positive base interval. */
Expr* interval_power_interval(const Expr* A, const Expr* B);

/* Companion set-operation builtins. */
Expr* builtin_intervalunion(Expr* res);
Expr* builtin_intervalintersection(Expr* res);
Expr* builtin_intervalmemberq(Expr* res);

#endif /* INTERVAL_H */

/* Mathilda — Rationalize.
 *
 * Converts inexact (Real / MPFR) numbers to nearby rationals, mirroring
 * Mathematica's Rationalize:
 *
 *   Rationalize[x]      → x converted to a nearby rational with small
 *                         denominator, but only if there is one within
 *                         |p/q − x| < 10^-4 / q^2; otherwise x is returned
 *                         unchanged.
 *   Rationalize[x, dx]  → the rational with the smallest denominator that
 *                         lies within dx of x (uses the standard "simplest
 *                         rational in interval" continued-fraction search).
 *   Rationalize[x, 0]   → forces conversion of any inexact x to rational
 *                         form, using a tolerance derived from the input's
 *                         precision (≈ ½·ulp).
 *
 * Threading:
 *   Rationalize maps over the internal structure of its argument:
 *     Rationalize[a + b]      ≡ Rationalize[a] + Rationalize[b]
 *     Rationalize[Complex[a,b]] ≡ Complex[Rationalize[a], Rationalize[b]]
 *   so e.g. Rationalize[1.2 + 6.7 x] becomes 6/5 + (67 x)/10.
 *
 * Two entry points are exposed: an `Expr*`-level builtin used by the
 * evaluator, and a lower-level `internal_rationalize_double` used by other
 * algorithms (Simplify, NSolve helpers, …) that already hold a `double`.
 */
#ifndef RATIONALIZE_H
#define RATIONALIZE_H

#include "expr.h"
#include <gmp.h>
#include <stdbool.h>

/* Mode selector for the core algorithm. */
typedef enum {
    /* No tolerance argument given: only return p/q if |p/q − x| < c/q^2
     * with c = 10^-4. Otherwise the caller should keep x unchanged. */
    RATIONALIZE_DEFAULT,
    /* Explicit positive tolerance dx > 0: return the rational with the
     * smallest denominator that lies in [x - dx, x + dx]. */
    RATIONALIZE_TOLERANCE,
    /* dx = 0: convert an inexact value to its precision-equivalent
     * smallest-denominator rational. Internally an ulp-derived tolerance is
     * used, matching Rationalize[x,0] in Mathematica. */
    RATIONALIZE_ZERO,
} RationalizeMode;

/* Core algorithm — purely numeric, no Expr involved.
 *
 * On entry: out_n / out_d may be uninitialised; on success they are
 * mpz_init'd and filled with the result in lowest terms.
 *
 * Returns true on success, false in RATIONALIZE_DEFAULT mode when no
 * convergent of x meets the c/q^2 bound. In the false case neither
 * out_n nor out_d is initialised, so the caller can leave them alone.
 *
 * This function is the "single-double" workhorse used both by the
 * Rationalize builtin and by other internal algorithms. */
bool internal_rationalize_double(double x, double dx, RationalizeMode mode,
                                 mpz_t out_n, mpz_t out_d);

/* Expression-level convenience: rationalize an entire expression by
 * threading over its structure. The caller-owned `e` is read but not
 * consumed. The returned Expr* is freshly allocated and owned by the
 * caller. Returns NULL only if `e` is NULL.
 *
 * Threading rules:
 *   - Inexact leaves (Real, MPFR) are converted via the core algorithm.
 *     In RATIONALIZE_DEFAULT mode, leaves whose value cannot meet the
 *     c/q^2 bound are returned unchanged.
 *   - Complex[re, im] is recursed into componentwise.
 *   - Integers, Rationals, Symbols are returned unchanged.
 *   - Other compound expressions: in TOLERANCE / ZERO modes, a fully
 *     numeric subexpression (per NumericQ) is first numericalised then
 *     rationalised — so `Rationalize[Sqrt[2], 0.01]` returns 99/70.
 *     Otherwise the function is rebuilt with each argument recursively
 *     rationalised, so `Rationalize[a + 1.5 x]` becomes `a + (3 x)/2`.
 *   - Hold / HoldForm / HoldComplete / HoldPattern / Unevaluated wrappers
 *     are passed through unchanged. */
Expr* internal_rationalize_expr(const Expr* e, double dx, RationalizeMode mode);

/* Builtin entry point. Implements Rationalize[x] and Rationalize[x, dx]. */
Expr* builtin_rationalize(Expr* res);

/* True iff the expression tree contains any inexact numeric leaf
 * (EXPR_REAL, or EXPR_MPFR when compiled with USE_MPFR). Walks the tree
 * including head and arguments of compound expressions. Used by the
 * exact-symbolic builtins (PolynomialGCD, Together, Apart, Factor, …) to
 * decide whether to force-rationalize their inputs. */
bool internal_contains_inexact(const Expr* e);

/* Convenience: true iff any argument of a function-call `res` contains an
 * inexact leaf. The head itself is not inspected. Returns false on a
 * non-FUNCTION input or one with no args. */
bool internal_args_contain_inexact(const Expr* res);

/* Force-rationalisation: convert every inexact leaf in `e` to a rational
 * so the downstream exact-arithmetic algorithms can run. Functionally
 * close to `Rationalize[expr, 0]` but with two pragmatic differences
 * tuned for internal use:
 *
 *   1. Per inexact leaf we try the convergent-bound algorithm first
 *      (Rationalize's no-tolerance form) — that recovers `1/9.` as
 *      `1/9` rather than the bit-exact `562949953421312/...`. Only
 *      when no convergent fits the bound do we fall back to the
 *      ½-ulp algorithm so the result is always a rational.
 *   2. Symbolic constants (Pi, E, ...) are NOT rationalised — they
 *      stay symbolic so e.g. `Together[Pi/x + 1.5]` keeps Pi exact and
 *      only converts the 1.5.
 *
 * The descent threads over compound expressions; Hold-forms are passed
 * through. Returned Expr is freshly allocated; `e` is read but not
 * consumed. */
Expr* internal_force_rationalize(const Expr* e);

/* Wrapper used by exact-symbolic builtins. If any argument of `res`
 * contains an inexact leaf, build a clone with each argument
 * force-rationalised and dispatch to `core(clone)`; on success the
 * returned expression is then numericalised so the caller observes
 * inexact-in / inexact-out semantics consistent with Mathematica.
 *
 * Idiomatic use at the top of a builtin:
 *
 *   if (internal_args_contain_inexact(res))
 *       return internal_rationalize_then_numericalize(res, builtin_xxx);
 *
 * Recursion bottoms out cleanly: the cloned `res` no longer contains
 * inexact leaves, so the inner `core(clone)` call falls through this
 * guard and runs the exact-arithmetic path.
 *
 * Ownership: matches the in-tree builtin convention — the evaluator
 * frees `res` after a non-NULL return, so this helper never touches
 * it. The internally-built clone is always freed by this helper. If
 * `core` returns NULL (cannot evaluate even the rationalised input),
 * the helper returns NULL too and the evaluator keeps the original
 * `res` unevaluated. */
Expr* internal_rationalize_then_numericalize(Expr* res,
                                             Expr* (*core)(Expr*));

/* Register the Rationalize symbol with the global symbol table. Called
 * from core_init. */
void rationalize_init(void);

#endif /* RATIONALIZE_H */

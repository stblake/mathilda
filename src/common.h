/* Mathilda — common preprocessing helpers.
 *
 * Several exact-symbolic builtins (Integrate, Solve, ...) need the same
 * "if the input has inexact (Real / MPFR) numbers, rationalise them so
 * the exact-arithmetic core can run, then numericalise the result so
 * the caller observes inexact-in / inexact-out semantics" preprocessing
 * step.  Rather than each routine re-implementing the detection +
 * rationalisation + post-numericalisation dance, this module exposes a
 * small shared API:
 *
 *   common_scan_inexact          - one pass; reports presence and the
 *                                  minimum precision (in bits) across
 *                                  every inexact leaf
 *   common_rationalize_input     - force-rationalise every inexact leaf
 *                                  at the given precision
 *   common_numericalize_result   - post-process an exact result back to
 *                                  that same precision
 *
 * The minimum-precision policy means a mixed input (say a machine-
 * precision Real `1.5` plus an MPFR value at 100 bits) propagates the
 * *lower* precision through to the answer: the answer is no more
 * precise than the least-precise input, which is the standard
 * inexact-arithmetic contract.  Pure-Real inputs get 53 bits (IEEE
 * double); pure-MPFR inputs get their explicit precision.
 *
 * The implementation delegates to `internal_force_rationalize_bits`
 * (rationalize.c) for leaf-level rationalisation and to `numericalize`
 * (numeric.c) for the post-processing.  Keeping the API here means
 * future callers (`NSolve`, `DSolve`, ...) plug in the same way.
 */
#ifndef COMMON_H
#define COMMON_H

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

/* True iff `e` is a function call whose head symbol equals `sym`.
 *
 * NULL-safe: returns false if `e` or `sym` is NULL or if `e` is not a
 * function with a symbol head.
 *
 * The comparison is *pointer equality* against the global symbol-name
 * intern table.  Every `EXPR_SYMBOL` carries a name pointer produced by
 * `intern_symbol`, so two heads are equal iff their pointers match.
 * Callers should pass one of the cached `SYM_*` constants from
 * `sym_names.h`; for runtime-supplied names, intern once with
 * `intern_symbol(name)` and pass the result.
 *
 * This is the single canonical head-check used throughout the system --
 * do not reintroduce strcmp-based head comparisons. */
bool head_is(const Expr* e, const char* sym);

/* Emit a Wolfram-style argument-count diagnostic on stderr for a builtin
 * `head` that was called with `got` arguments but expects an argument count
 * in the inclusive range [min, max].
 *
 *   - min == max            -> fixed arity   ("Sin::argx" / "Head::argrx")
 *   - max == SIZE_MAX       -> min-or-more   ("Head::argm")
 *   - max == min + 1        -> two choices   ("Head::argt": "1 or 2 ...")
 *   - otherwise             -> a range       ("Head::argb": "between m and n")
 *
 * The message text mirrors Mathematica exactly, e.g.
 *   Sin::argx: Sin called with 3 arguments; 1 argument is expected.
 *
 * Always returns NULL so an arity guard reads as a single statement:
 *   if (argc != 1)
 *       return builtin_arg_error("Sin", argc, 1, 1);
 * Leaving the call unevaluated (NULL) matches Mathematica, which prints the
 * message and returns the expression unchanged. */
Expr* builtin_arg_error(const char* head, size_t got, size_t min, size_t max);

/* Result of common_scan_inexact.
 *
 *   has_inexact: true iff `e` contains any inexact numeric leaf
 *                (EXPR_REAL, or EXPR_MPFR when compiled with USE_MPFR).
 *   min_bits:    minimum bit-precision across all inexact leaves.
 *                Real leaves contribute 53 (IEEE 754 double);
 *                MPFR leaves contribute their mpfr_get_prec().
 *                Undefined (zero) when has_inexact is false. */
typedef struct {
    bool has_inexact;
    long min_bits;
} CommonInexactInfo;

/* One-pass scan: walks heads and arguments of compound expressions,
 * collecting both the presence flag and the minimum bit-precision
 * across all inexact leaves.  NULL-safe (returns {false, 0}). */
CommonInexactInfo common_scan_inexact(const Expr* e);

/* Return a freshly allocated copy of `e` with every inexact numeric
 * leaf replaced by an exact rational at the given precision.  `bits`
 * is the precision in bits (typically the `min_bits` reported by
 * common_scan_inexact); the rationalisation uses a tolerance of
 * 2^(-bits), so the result is no more precise than the input.
 * Symbolic constants (Pi, E, ...) and Hold-forms pass through
 * unchanged.  Caller owns the returned Expr*. */
Expr* common_rationalize_input(const Expr* e, long bits);

/* Round-trip an exact result back to floating-point at the given
 * precision.  Used as the inverse of common_rationalize_input after a
 * routine such as Solve / Integrate has produced an exact symbolic
 * answer from a rationalised input, so the user sees inexact-in /
 * inexact-out semantics.
 *
 *   bits <= 53        →  IEEE 754 double (EXPR_REAL).
 *   bits >  53        →  MPFR at exactly `bits` bits (when USE_MPFR is
 *                        defined; falls back to machine otherwise).
 *
 * The traversal is structural: lists, rules, and other compound
 * expressions are recursed into so e.g. `{{x -> 1/2}}` numericalises
 * to `{{x -> 0.5}}` while leaving unknown symbols alone.
 *
 * Caller owns the returned Expr*. */
Expr* common_numericalize_result(const Expr* e, long bits);

/* ---------------------------------------------------------------------- *
 * Per-method head aliases                                                 *
 *                                                                         *
 * A builtin that runs a cascade of strategies behind `Method -> m` also   *
 * exposes each setting as a head of its own -- Limit`Series[f, x -> a] is *
 * exactly Limit[f, x -> a, Method -> "Series"], NLimit`EulerSum[...] is   *
 * NLimit[..., Method -> "EulerSum"], and so on.  This is the shared       *
 * plumbing behind those heads.                                            *
 *                                                                         *
 * `res` is rebuilt as head[args..., Method -> method] and handed to       *
 * `impl`, the standard builtin for `head`.  Any Method option the caller  *
 * supplied among the trailing options (everything past `n_positional`) is *
 * dropped: the head already names the method, and honouring a second,     *
 * possibly conflicting, setting would make                                *
 * Limit`Series[f, s, Method -> "Gruntz"] ambiguous.  Every other option   *
 * (Direction, WorkingPrecision, ...) is forwarded untouched.              *
 *                                                                         *
 * `impl` is called directly rather than through the evaluator, so the     *
 * builtin ownership contract puts this function in the evaluator's seat:  *
 * the constructed call is freed here, and the caller receives either a    *
 * fresh result it owns or NULL (leaving the alias call unevaluated, which *
 * is the same "this method does not apply" answer the Method option       *
 * gives).  `res` itself is never freed.                                   *
 * ---------------------------------------------------------------------- */
Expr* common_method_alias(Expr* res, const char* head, size_t n_positional,
                          const char* method, Expr* (*impl)(Expr*));

#endif /* COMMON_H */

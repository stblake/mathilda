/*
 * Mathilda — PossibleZeroQ / hybrid symbolic-numeric zero recognition.
 *
 * A heuristic decision procedure for "is this expression equivalent to 0?".
 * The general problem is undecidable (Richardson 1968), so PossibleZeroQ
 * combines several decidable / probabilistic sub-procedures:
 *
 *   Stage 0  — Structural shortcuts (literal 0, Complex[0,0], …).
 *   Stage 1  — Rational normalization (Together / Cancel / Expand +
 *              is_zero_poly). Decides every identity in Q(x_1,…,x_n).
 *   Stage 2  — Numeric precision ladder on closed-form numeric inputs.
 *              Compares |numericalize(e, p)| against a cancellation-scaled
 *              ambiguity threshold; bumps precision (machine → 200 → 500
 *              → 1000 bits) to distinguish true zero from catastrophic
 *              cancellation.
 *   Stage 3  — Schwartz–Zippel probabilistic identity test: substitute
 *              free symbols with random rationals (drawn from Q[i]) and
 *              recurse into Stage 2. Requires k=4 independent confirmations.
 *
 * See ZERO_RECOGNISE_PLAN.md for the full design rationale and literature.
 */
#ifndef MATHILDA_ZERO_TEST_H
#define MATHILDA_ZERO_TEST_H

#include "expr.h"

typedef enum {
    ZERO_TEST_FALSE   = 0,  /* proved (or strongly believed) non-zero       */
    ZERO_TEST_TRUE    = 1,  /* proved zero (rational) or numerically zero   */
    ZERO_TEST_UNKNOWN = 2   /* heuristic exhausted; caller decides policy   */
} ZeroTestResult;

/* Three-valued internal decision procedure. The input `e` is read-only and
 * is NOT consumed. Callers (Equal, Simplify, integration) that need to
 * distinguish "unknown" from "true zero" should use this entry. */
ZeroTestResult zero_test_decide(const Expr* e);

/* PossibleZeroQ[expr] builtin. Follows the standard Mathilda contract:
 *   - Consumes ownership of `res` indirectly (evaluator frees it).
 *   - Returns a fresh True/False symbol, or NULL on arity mismatch.
 *
 * Collapses ZERO_TEST_UNKNOWN to True (matching Mathematica's
 * PossibleZeroQ::ztest1 "assume zero when uncertain" behaviour). */
Expr* builtin_possible_zero_q(Expr* res);

/* Registers PossibleZeroQ in the symbol table. Called from core_init(). */
void zero_test_init(void);

#endif /* MATHILDA_ZERO_TEST_H */

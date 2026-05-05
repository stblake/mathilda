#ifndef TRIG_CANON_H
#define TRIG_CANON_H

#include "expr.h"
#include <stddef.h>

/* A (base, exponent) pair as accumulated by builtin_times after grouping. */
typedef struct {
    Expr* base;
    Expr* exponent;
} BasePower;

/* Canonicalize trig and hyperbolic ratios across an array of (base, exponent)
 * groups, in place.
 *
 * Each contributing factor is mapped to a (Sin^a * Cos^b) decomposition, and
 * factors that share the same argument are combined and re-emitted using the
 * shortest naming:
 *
 *     Sin[x] * Cos[x]^-1   -> Tan[x]
 *     Cos[x]^-1            -> Sec[x]
 *     Sin[x] * Csc[x]      -> 1   (group removed)
 *     Tan[x] * Cos[x]      -> Sin[x]
 *     Sec[x] * Tan[x]      -> Sec[x] * Tan[x]   (no shorter form)
 *     Sin[x]^2 / Cos[x]^2  -> Tan[x]^2
 *
 * The same set of rewrites is applied to the hyperbolic family
 * {Sinh, Cosh, Tanh, Coth, Sech, Csch}, independently from the trig family.
 *
 * Only fires when the exponent of a candidate group is a plain EXPR_INTEGER.
 * Non-integer or symbolic exponents are left alone. The total group_count
 * never grows under this transform: solo rewrites preserve count, and
 * multi-contributor merges always emit at most as many entries as they consume
 * (each bucket emits 0, 1, or 2 entries).
 */
void trig_canon_groups(BasePower* groups, size_t* group_count);

/* Standalone Power-level reciprocal rewrite. The parser collapses things like
 * `1/Cos[x]` into `Power[Cos[x], -1]` without going through Times, so the
 * Times-level grouping pass never sees them. This helper handles that case.
 *
 *     Power[Cos[x], -1]  -> Sec[x]
 *     Power[Sin[x], -2]  -> Csc[x]^2
 *     Power[Tan[x], -1]  -> Cot[x]
 *     Power[Sech[x], -3] -> Cosh[x]^3
 *
 * Returns a freshly-allocated Expr on a successful rewrite, or NULL when the
 * base is not a trig/hyp ratio function, the exponent is not a negative
 * integer, or no shorter form exists (e.g., positive exponents are already
 * canonical and left alone).
 */
Expr* trig_canon_power(Expr* base, int64_t exp);

/* Reentrant suppress counter. Increment to disable trig_canon's rewrites
 * (both the Times-level pass and the Power-level reciprocal rewrite) within
 * the dynamic extent of an operation, then decrement to restore.
 *
 * Why this exists: trigsimp's pipeline (TrigReduce, TrigFactor, ExpToTrig,
 * Simplify) internally rewrites Tan -> Sin/Cos, performs algebraic work on
 * the Sin/Cos polynomial form, then converts back. Without suppression, the
 * very first step is undone immediately by the Times-level canonicalizer
 * (Sin/Cos collapses straight back to Tan), so the pipeline cannot make
 * progress. The suppress region lets the pipeline operate on the explicit
 * Sin/Cos form internally; once it returns, the outer evaluator re-evaluates
 * the result and the canonicalizer fires normally on the final answer.
 *
 * Pairs MUST be balanced. Calls nest: inc/inc/dec/dec stays suppressed until
 * the outermost dec brings the counter back to zero. */
void trig_canon_suppress_inc(void);
void trig_canon_suppress_dec(void);

#endif

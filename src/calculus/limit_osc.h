/*
 * limit_osc.h -- Oscillatory normal form for Limit[f, x -> +/-Infinity].
 *
 * A decision layer for limits of *amplitude-modulated oscillations*: the
 * shapes where a trigonometric factor whose argument diverges is multiplied
 * by an algebraic / exp-log amplitude. The layer either
 *
 *   - returns the limit (when the oscillation provably washes out, or is
 *     provably dominated by a non-oscillatory part), or
 *   - returns Indeterminate (when the oscillation provably survives, so no
 *     limit exists), or
 *   - returns NULL (abstain -- the caller keeps walking its cascade).
 *
 * All three verdicts are theorems; the module never guesses. See limit_osc.c
 * for the normal form, the hypotheses each rule checks, and the proofs.
 *
 * The module is deliberately decoupled from limit.c's private LimitCtx: sub-
 * limits (of amplitudes and of amplitude ratios) are requested through a
 * caller-supplied callback so the caller keeps its own recursion-depth and
 * Direction bookkeeping.
 */

#ifndef MATHILDA_LIMIT_OSC_H
#define MATHILDA_LIMIT_OSC_H

#include "expr.h"

/* Compute Limit[g, x -> point] in the caller's context.
 *
 * `g` is borrowed (the callback must not free or mutate it). The return is
 * either a freshly-allocated Expr* owned by the caller of the callback, or
 * NULL when the sub-limit could not be determined. */
typedef Expr* (*LimitOscSubFn)(Expr* g, void* ctx);

/* Decide Limit[f, x -> point] for the oscillatory class.
 *
 * `point` must be Infinity or -Infinity; anything else abstains immediately.
 * A finite point is the caller's job: substitute x = a +/- 1/t and ask for
 * t -> +Infinity (limit.c's layer_oscillatory does exactly that).
 * `f`, `x` and `point` are borrowed. Returns a freshly-allocated Expr* (the
 * limit value, or the symbol Indeterminate when no limit exists) or NULL to
 * abstain. */
Expr* limit_oscillatory(Expr* f, Expr* x, Expr* point,
                        LimitOscSubFn sub, void* subctx);

#endif /* MATHILDA_LIMIT_OSC_H */

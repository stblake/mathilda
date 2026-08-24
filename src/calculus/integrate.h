/* integrate.h
 *
 * Public entry points for the rational-function integrator.  The
 * top-level dispatcher `Integrate[f, x]` lives in the System` context
 * and routes to `Integrate`BronsteinRational[f, x]` (implemented in
 * intrat.c) when its input is a polynomial or a rational function in x.
 *
 * Phase 1 of the BronsteinRational port (see plans/INTEGRATE_PLAN.md).  The
 * non-rational fallback is the identity (returns the call unevaluated)
 * so this file remains a tiny shim — all real work is in intrat.c.
 */

#ifndef MATHILDA_INTEGRATE_H
#define MATHILDA_INTEGRATE_H

#include "expr.h"

/* `Integrate[f, x]` — System` builtin.  Validates that `f` is a
 * polynomial or a rational function in `x`; on success forwards to
 * `Integrate`BronsteinRational[f, x]`.  Returns NULL (leaving the call
 * unevaluated) for any other input. */
Expr* builtin_integrate(Expr* res);

/* Nesting depth of the public Integrate builtin's method cascade.  0 outside
 * any integration; 1 while the outermost (user-facing) Integrate[f, x] runs its
 * cascade; >= 2 inside a sub-integral spawned by an internal substitution
 * (DerivativeDivides' u = Log[x], etc.) that re-enters Integrate via the
 * evaluator.  The deep RischTranscendental decision stage reads this to fire the
 * user-facing Integrate::nonelem diagnostic only for the ORIGINAL integrand
 * (depth <= 1), never for an internal recursion variable. */
extern int g_integrate_depth;

/* Registers `Integrate` in the symbol table along with its docstring
 * and attributes.  Also calls `intrat_init()` to register every
 * `Integrate`...` package builtin so they are available before the
 * REPL accepts user input. */
void integrate_init(void);

#endif /* MATHILDA_INTEGRATE_H */

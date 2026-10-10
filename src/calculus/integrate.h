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

/* Suppression counter for the user-facing Integrate::nonelem diagnostic.  When
 * > 0, the diagnostic is not printed even at the outermost frame.  A caller that
 * integrates *speculatively* — where a non-elementary result means "this method
 * declines", not "warn the user about their input" — brackets its ds_integrate
 * calls with g_integrate_quiet++/-- so those declines stay silent.  The Kovacic
 * solver (recovery weight, second solution, apparent-singularity integrals) is
 * the first such caller: a RootSum-valued log part there is a decline, not a
 * report about the ODE the user actually typed. */
extern int g_integrate_quiet;

/* Suppression counter for the heavy TAIL stages of the Automatic cascade
 * (ParallelMixedTower, GoursatAlgebraic, ParallelMixedSpecial).  When > 0, the
 * cascade runs only its cheap, elementary stages and otherwise declines rather
 * than paying the seconds-long tower / special-function search.  A stage that
 * integrates a SUB-integrand speculatively — where "did not close cheaply" is a
 * decline, not a reason to grind — brackets its recursive Integrate calls with
 * g_integrate_no_special++/--, so a non-closing sub-integral returns promptly and
 * the stage declines fast, letting the OUTER cascade (which the flag does not
 * touch, being back to 0 there) take its own turn at the tail stages.  The
 * LogByParts stage is the first such caller: its V = INT K and W = INT V*(g'/g)
 * reductions close at a cheap stage on the Log*trig class it targets, and a case
 * where they do not must not cost two passes through ParallelMixedSpecial. */
extern int g_integrate_no_special;

/* Emits the user-facing `Integrate::nonelem` diagnostic for the original
 * integrand `f` (in variable `x`), shared by every method that can PROVE the
 * integrand has no elementary antiderivative (RischTranscendental via its field
 * decision; ParallelMixedTower via a {"not elementary", ...} certificate).  The
 * caller must have already established the proof; this routine handles only the
 * gating (silent for internal recursion / speculative integration, and at most
 * once per top-level cascade) and the message.  Borrows `f` and `x`. */
void integrate_announce_nonelementary(Expr* f, Expr* x);

/* Registers `Integrate` in the symbol table along with its docstring
 * and attributes.  Also calls `intrat_init()` to register every
 * `Integrate`...` package builtin so they are available before the
 * REPL accepts user input. */
void integrate_init(void);

#endif /* MATHILDA_INTEGRATE_H */

/* integrate_risch_macsyma.h — Risch integrator ported from Maxima.
 *
 * `Integrate`RischMacsyma[f, x]` is a faithful port of the algorithm
 * structure of Maxima's src/risch.lisp: a recursive Risch integrator
 * organised around a differential transcendental tower over a single
 * integration variable, with dedicated rational, logarithmic,
 * exponential, and (flag-gated) special-function cases.
 *
 * Unlike Maxima, arithmetic is grounded in Mathilda's existing
 * Expr/poly/rat machinery (Together, Apart, Factor, PolynomialGCD,
 * RowReduce, D) rather than a reimplementation of Maxima's CRE
 * polynomial-list representation.  The rational case delegates to
 * `Integrate`BronsteinRational`.
 *
 * It is reachable via the explicit backtick form and via
 * `Integrate[f, x, Method -> "RischMacsyma"]`, and it is inserted into
 * the Integrate Automatic cascade after `Integrate`RischNorman`.  As a
 * decision procedure, every branch is correct by construction: each case
 * fires only behind an exact structural certificate that already proves
 * the closed form it emits, so the result is NOT checked by
 * differentiation (Risch integration is not a guess-and-verify search).
 *
 * This is deliberately distinct from `Integrate`RischNorman`, which is the
 * *parallel* Risch (Norman/pmint) heuristic; RischMacsyma is the recursive
 * Risch algorithm and never falls back on the parallel-Risch engine.
 *
 * Memory contract follows the standard Mathilda BuiltinFunc rule: the
 * caller (evaluator) owns `res`; on success the builtin returns a
 * freshly-allocated Expr*, on failure NULL (the call bubbles symbolic).
 *
 * Implementation phases (see the approved plan):
 *   Phase 1: scaffold + dispatch + rational case + verification gate.
 *   Phase 2: differential tower + logarithmic case.
 *   Phase 3: exponential case (Risch differential equation).
 *   Phase 4: shared linear solver + polynomial-part recursion + caps.
 *   Phase 5: special functions (Erf/PolyLog, gated) + trig front-end.
 */

#ifndef MATHILDA_INTEGRATE_RISCH_MACSYMA_H
#define MATHILDA_INTEGRATE_RISCH_MACSYMA_H

#include "expr.h"

/* Public package entry: corresponds to Maxima's $risch(exp, var).
 * Returns a freshly-allocated antiderivative Expr* on success, or NULL
 * on failure (the call bubbles back as unevaluated). */
Expr* builtin_rischmacsyma(Expr* res);

/* Register the Integrate`RischMacsyma package symbol in the global
 * symbol table.  Called from integrate_init() during core_init().
 * Idempotent. */
void integrate_risch_macsyma_init(void);

#endif /* MATHILDA_INTEGRATE_RISCH_MACSYMA_H */

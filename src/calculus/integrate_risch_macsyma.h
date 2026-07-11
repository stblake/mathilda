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

#include <stdbool.h>
#include "expr.h"

/* Bronstein RdeBoundDegree leading-degree bound for the Risch differential equation
 * D[q] + f q = p in a single monomial v.  Returns an upper bound on deg_v(q).
 *   dpv          deg_v(p)  (as deg(numerator) - deg(denominator))
 *   dfv          deg_v(f)  (likewise)
 *   deriv_lowers true for a deriv-lowering monomial (base x under d/dx, or a logarithmic
 *                kernel); false for a deriv-preserving one (an exponential kernel).
 *   m_res        the Bronstein resonance integer for the leading-coefficient
 *                cancellation sub-case (nonnegative), or -1 when none applies; the bound
 *                is widened monotonically to max(naive, m_res) in the cancellation config.
 * See the definition in integrate_risch_macsyma.c for the full derivation and a proof
 * that both cancellation configurations are pre-empted (so m_res == -1) on every
 * reachable call.  Exposed so the degree arithmetic — including the resonance widening —
 * can be unit-tested directly without constructing an integrand. */
long rm_rde_var_bound(long dpv, long dfv, bool deriv_lowers, long m_res);

/* Public package entry: corresponds to Maxima's $risch(exp, var).
 * Returns a freshly-allocated antiderivative Expr* on success, or NULL
 * on failure (the call bubbles back as unevaluated). */
Expr* builtin_rischmacsyma(Expr* res);

/* Register the Integrate`RischMacsyma package symbol in the global
 * symbol table.  Called from integrate_init() during core_init().
 * Idempotent. */
void integrate_risch_macsyma_init(void);

#endif /* MATHILDA_INTEGRATE_RISCH_MACSYMA_H */

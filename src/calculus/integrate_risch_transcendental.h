/* integrate_risch_transcendental.h — recursive transcendental Risch integrator.
 *
 * `Integrate`RischTranscendental[f, x]` is the recursive Risch decision
 * procedure for transcendental elementary functions (the Bronstein/Roach
 * lineage): a recursive integrator organised around a differential
 * transcendental tower over a single integration variable, with dedicated
 * rational, logarithmic, exponential, and (flag-gated) special-function
 * cases.
 *
 * Arithmetic is grounded in Mathilda's existing Expr/poly/rat machinery
 * (Together, Apart, Factor, PolynomialGCD, RowReduce, D) — each kernel is
 * aliased to a fresh polynomial variable and manipulated with those
 * primitives; there is no separate polynomial-list representation.  The
 * rational base case delegates to `Integrate`BronsteinRational`.
 *
 * It is reachable via the explicit backtick form and via
 * `Integrate[f, x, Method -> "RischTranscendental"]`, and it is inserted into
 * the Integrate Automatic cascade after `Integrate`RischNorman`.  As a
 * decision procedure, every branch is correct by construction: each case
 * fires only behind an exact structural certificate that already proves
 * the closed form it emits, so the result is NOT checked by
 * differentiation (Risch integration is not a guess-and-verify search).
 *
 * This is deliberately distinct from `Integrate`RischNorman`, which is the
 * *parallel* Risch (Norman/pmint) heuristic; RischTranscendental is the recursive
 * Risch algorithm and never falls back on the parallel-Risch engine.
 *
 * Memory contract follows the standard Mathilda BuiltinFunc rule: the
 * caller (evaluator) owns `res`; on success the builtin returns a
 * freshly-allocated Expr*, on failure NULL (the call bubbles symbolic).
 */

#ifndef MATHILDA_INTEGRATE_RISCH_TRANSCENDENTAL_H
#define MATHILDA_INTEGRATE_RISCH_TRANSCENDENTAL_H

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
 * See the definition in integrate_risch_transcendental.c for the full derivation and a proof
 * that both cancellation configurations are pre-empted (so m_res == -1) on every
 * reachable call.  Exposed so the degree arithmetic — including the resonance widening —
 * can be unit-tested directly without constructing an integrand. */
long rt_rde_var_bound(long dpv, long dfv, bool deriv_lowers, long m_res);

/* Public package entry for `Integrate`RischTranscendental[f, x]`.
 * Returns a freshly-allocated antiderivative Expr* on success, or NULL
 * on failure (the call bubbles back as unevaluated). */
Expr* builtin_rischtranscendental(Expr* res);

/* Enable/disable the constant-base debasing pre-pass (a^u -> E^(u Log a)).
 * Returns the PREVIOUS setting so callers can save/restore.  The Integrate
 * inexact path disables it around the rationalised cascade: a float base
 * (2.71828^x) is rationalised to a large exact fraction upstream, and debasing
 * that fraction would manufacture a spuriously "exact" Log-based closed form
 * that then numericalises to a degraded exponent — inexact inputs must instead
 * observe the conservative inexact-in / inexact-out contract. */
bool rt_transcendental_set_debase(bool enabled);

/* Register the Integrate`RischTranscendental package symbol in the global
 * symbol table.  Called from integrate_init() during core_init().
 * Idempotent. */
void integrate_risch_transcendental_init(void);

#endif /* MATHILDA_INTEGRATE_RISCH_TRANSCENDENTAL_H */

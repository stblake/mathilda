/* risch_hypertangent.h — the hypertangent case (Bronstein §5.10).
 *
 * A hypertangent monomial t over k satisfies Dt/(t^2+1) = a in k (t = tan(∫a)):
 * a NONLINEAR monomial whose only special irreducible is t^2+1.  Bronstein §5.10
 * integrates over such an extension DIRECTLY, keeping tan/arctan real, WITHOUT
 * rewriting to complex exponentials.
 *
 * This module implements the polynomial-part algorithm IntegrateHypertangent-
 * Polynomial (§5.10, p.167): given p in k[t] it returns q in k[t] and c in k with
 *
 *     p - D[q] - c D(t^2+1)/(t^2+1)  in k                                (*)
 *
 * so that ∫ p = q + c log(t^2+1) + ∫(element of k).  Since D(t^2+1)/(t^2+1) =
 * 2 a t, this is the tangent analogue of the exponential/logarithmic polynomial
 * cases.  If Dc != 0 the polynomial p has no elementary integral over k
 * (the log-term coefficient must be a constant, by Liouville's theorem).
 *
 * Exposed as an internal builtin:
 *   Risch`IntegrateHypertangentPolynomial[p, t, deriv] -> {q, c}
 * where deriv = {x -> 1, ..., t -> Dt, ...} with Dt = a (t^2+1).
 */

#ifndef MATHILDA_RISCH_HYPERTANGENT_H
#define MATHILDA_RISCH_HYPERTANGENT_H

#include <stdbool.h>
#include "expr.h"
#include "risch_field.h"

/* IntegrateHypertangentPolynomial (§5.10).  For a hypertangent monomial t of the
 * derivation d and p in k[t], writes owned q in k[t] and c in k satisfying (*).
 * Returns false if t is not a hypertangent monomial of d. */
bool risch_integrate_hypertangent_poly(const Expr* p, const Expr* t,
                                       const RischDeriv* d, Expr** q, Expr** c);

/* ResidueReduce (§5.6, the residue criterion).  For the simple part h in k(t),
 * writes owned g2 (a sum of logs) and *beta: true with h - D[g2] in k[t] when
 * every Rothstein-Trager residue is a constant of d, or false (g2 = 0) when a
 * residue is non-constant (h has no elementary integral over k(t)). */
bool risch_residue_reduce(const Expr* h, const Expr* t, const RischDeriv* d,
                          Expr** g2, bool* beta);

/* IntegrateHypertangent (§5.10, p.172 — the full driver).  For a hypertangent
 * monomial t of the derivation given by the rules List `deriv` and f in k(t),
 * writes owned g and *beta: true with f - D[g] in k, or false when f - D[g] has
 * no elementary integral over k(t).  Composes HermiteReduce, ResidueReduce, and
 * the reduced and polynomial hypertangent cases. */
bool risch_integrate_hypertangent(const Expr* f, const Expr* t,
                                  const Expr* deriv, Expr** g, bool* beta);

/* IntegrateHypertanh (the hyperbolic tangent case, special t^2-1).  For a
 * hyperbolic tangent monomial t (Dt = a(t^2-1)) of the derivation `deriv` and f
 * in k(t), writes owned g and *beta: true with f - D[g] in k, or false when
 * f - D[g] has no elementary integral over k(t).  The special t^2-1 splits over
 * k, so the reduced case decouples into two real base Risch DEs (no C(i)). */
bool risch_integrate_hypertanh(const Expr* f, const Expr* t,
                               const Expr* deriv, Expr** g, bool* beta);

/* Register the hypertangent builtins.  Called from integrate_init(). */
void integrate_risch_hypertangent_init(void);

#endif /* MATHILDA_RISCH_HYPERTANGENT_H */

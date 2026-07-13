/* risch_hermite.h — the Hermite reduction over a monomial extension.
 *
 * Bronstein, *Symbolic Integration I*, 2nd ed., §5.3 (Theorem 5.3.1, the
 * quadratic HermiteReduce algorithm, p.139).
 *
 * Given a derivation D on k[t] and f in k(t), the Hermite reduction computes
 * g, h, r in k(t) with
 *
 *     f = D[g] + h + r,   h SIMPLE (squarefree normal denominator),
 *                          r REDUCED (in k<t>: no normal poles left).
 *
 * It rewrites the normal part of f (the part with normal poles) as a derivative
 * D[g] plus a simple remainder, using only the extended Euclidean algorithm
 * (Diophantine solves) in k[t] — no factorization beyond squarefree.  The
 * simple part h then flows to the residue criterion / Rothstein-Trager log
 * part, and the reduced part r to the polynomial reduction.
 *
 * This is the literal Bronstein algorithm over the tower derivation, replacing
 * the earlier undetermined-coefficient Hermite ANSATZ (which restricted the
 * numerator coefficients to polynomials of a heuristic x-degree).  It is exact
 * and works with arbitrary rational k = C(x, t_1, ..., t_{i-1}) coefficients.
 *
 * Exposed as an internal builtin:
 *   Risch`HermiteReduce[f, t, deriv] -> {g, h, r}
 * where deriv = List[Rule[var, Dvar], ...] names the base variable (Dvar = 1)
 * and each monomial variable with its derivative.
 */

#ifndef MATHILDA_RISCH_HERMITE_H
#define MATHILDA_RISCH_HERMITE_H

#include <stdbool.h>
#include "expr.h"
#include "risch_field.h"

/* HermiteReduce (§5.3).  Writes owned g, h, r in k(t) with f = D[g] + h + r,
 * h simple and r reduced.  Returns false only on malformed input (t not a
 * monomial variable of d). */
bool risch_hermite_reduce(const Expr* f, const Expr* t, const RischDeriv* d,
                          Expr** g, Expr** h, Expr** r);

/* Register the Hermite-reduction builtin.  Called from integrate_init(). */
void integrate_risch_hermite_init(void);

#endif /* MATHILDA_RISCH_HERMITE_H */

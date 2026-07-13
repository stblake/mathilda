/* risch_coupled.h — coupled differential systems + hypertangent reduced case.
 *
 * Bronstein, *Symbolic Integration I*, 2nd ed., Chapter 8 (coupled differential
 * systems) and §5.10 (the hypertangent case, reduced part).
 *
 * The coupled 2x2 real system over a differential field (k, D) (Bronstein
 * Eq. 8.1, the a = -1 tangent case):
 *
 *     D y1 + f1 y1 - f2 y2 = g1
 *     D y2 + f2 y1 + f1 y2 = g2                                (f1,f2,g1,g2 in k)
 *
 * is equivalent, via y = y1 + y2 i, f = f1 + f2 i, g = g1 + g2 i, to the SINGLE
 * Risch differential equation  D y + f y = g  over k(i) = k(sqrt(-1))
 * (Bronstein §8.1, "the trivial route").  For the base field k = C(x) the
 * solution components y1, y2 come back real (in C(x)), so this reduction is
 * both correct and real-valued — it is what IntegrateHypertangentReduced needs.
 * The full Chapter-8 real recursion (CoupledDECancelTan), which avoids k(i) when
 * k itself carries further tangent monomials, is a later refinement.
 *
 * IntegrateHypertangentReduced (§5.10, p.169): for a hypertangent monomial t
 * (Dt/(t^2+1) = eta in k) and p in k(t) whose only poles are at t^2+1, peels the
 * poles one multiplicity at a time by solving the coupled system Eq. (5.20)
 *     (Dc; Dd) + [[0, -2 m eta],[2 m eta, 0]] (c; d) = (a; b),
 * returning q in k(t) with p - D[q] in k[t], or proving non-elementarity.
 *
 * Exposed as internal builtins:
 *   Risch`CoupledDESystem[f1, f2, g1, g2, x]         -> {y1, y2} | unevaluated
 *   Risch`IntegrateHypertangentReduced[p, t, deriv]  -> {q, beta}
 */

#ifndef MATHILDA_RISCH_COUPLED_H
#define MATHILDA_RISCH_COUPLED_H

#include <stdbool.h>
#include "expr.h"
#include "risch_field.h"

/* CoupledDESystem (a = -1), over the base field C(x).  Solves the 2x2 system
 * above for owned y1, y2 in C(x).  Returns false ("no solution") when the
 * equivalent Risch DE over C(i)(x) has no rational solution. */
bool risch_coupled_desystem(const Expr* f1, const Expr* f2,
                            const Expr* g1, const Expr* g2, const Expr* x,
                            Expr** y1, Expr** y2);

/* IntegrateHypertangentReduced (§5.10).  For a hypertangent monomial t of the
 * derivation d and p in k(t) with poles only at t^2+1, writes owned q in k(t)
 * with p - D[q] in k[t]; *beta is set to true, or to false when p has no
 * elementary integral over k(t) (a coupled system had no solution).  Returns
 * false if t is not a hypertangent monomial of d, or the input is malformed. */
bool risch_integrate_hypertangent_reduced(const Expr* p, const Expr* t,
                                          const RischDeriv* d, Expr** q, bool* beta);

/* Register the coupled-system / hypertangent-reduced builtins. */
void integrate_risch_coupled_init(void);

#endif /* MATHILDA_RISCH_COUPLED_H */

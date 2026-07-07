/*
 * integrate_residue.h -- Definite integration by the residue theorem.
 *
 * A correct-by-construction, NIntegrate-gated method for the classical families
 * of improper and periodic real integrals that complex analysis dispatches by
 * summing residues over the poles enclosed by a standard contour:
 *
 *   Family A -- rational integrands on (-Inf, Inf):
 *        Integrate[1/(1+x^4), {x, -Infinity, Infinity}]  ->  Pi/Sqrt[2]
 *        (value = 2 Pi I * Sum of residues in the upper half-plane).
 *
 *   Family B -- Fourier / Jordan integrands R(x) * {Cos, Sin, Exp[I .]}(a x) on
 *        (-Inf, Inf):
 *        Integrate[Cos[x]/(1+x^2), {x, -Infinity, Infinity}]  ->  Pi/E
 *        (J = 2 Pi I * Sum_UHP Res[R Exp[I a x]]; Cos -> Re[J], Sin -> Im[J]).
 *
 *   Family C -- rational-in-{Sin,Cos} integrands over a full period (0,2Pi) or
 *        (-Pi,Pi), via z = Exp[I x] on the unit circle:
 *        Integrate[1/(2+Cos[x]), {x, 0, 2 Pi}]  ->  2 Pi/Sqrt[3]
 *        (value = 2 Pi I * Sum of residues inside the unit disk).
 *
 *   Principal value -- a simple pole on the real axis (families A/B) contributes
 *        a HALF residue (Pi I Res), enabling e.g.
 *        Integrate[Sin[x]/x, {x, -Infinity, Infinity}]  ->  Pi.
 *
 *   Half-line [0, Inf) -- when f is even, Integrate[f, {x, 0, Infinity}] is half
 *        the full-line value: Integrate[1/(1+x^4), {x, 0, Infinity}] -> Pi/(2 Sqrt[2]).
 *
 * Each family's value is correct BY CONSTRUCTION once its structural gates hold,
 * so there is no numeric quadrature crosscheck; a result that does not close to a
 * scalar (surviving Root/radical), or whose imaginary part does not vanish for a
 * real-valued family (a mis-classified pole), returns NULL so the existing
 * Newton-Leibniz path takes over.  Reached via Integrate[f, {x, a, b}]
 * (Automatic, before Newton-Leibniz) or Integrate[f, {x, a, b},
 * Method -> "Residue"], and directly via Integrate`ContourResidue[f, {x, a, b}].
 */

#ifndef MATHILDA_INTEGRATE_RESIDUE_H
#define MATHILDA_INTEGRATE_RESIDUE_H

#include <stdbool.h>
#include "expr.h"

/*
 * Master entry for the Integrate dispatcher's definite real-spec path.  f, x, a,
 * b are borrowed (not consumed); x must be a symbol.  Returns a freshly-allocated
 * Expr* (the closed-form value) when a recognizer fires and its value closes to a
 * scalar, or NULL to leave the integral for the Newton-Leibniz path (no family
 * matched, or closure failed).
 *
 * `diverges` (may be NULL) is an out-flag, set to true only when a family
 * CONCLUSIVELY determines the ordinary integral does not converge -- a genuine
 * (non-removable) pole on the integration contour, i.e. the real axis for the
 * rational family or the unit circle for the trig family.  The dispatcher uses
 * it to emit Integrate::idiv and stop, rather than silently falling through.  It
 * is left false for a merely-unrecognized integrand (an undecidable/symbolic
 * pole, closure failure, etc.), which stays eligible for the Newton-Leibniz path.
 *
 * `assumptions` (may be NULL) is the borrowed value of the caller's
 * `Assumptions -> ...` option -- a fact expression (And / List of relations such
 * as `a > 0`, `0 < a < 1`, `n > 1`) constraining the integrand's free
 * parameters.  When present, the recognizers classify poles and frequencies that
 * depend on symbolic parameters by reading their sign at a single generic point
 * of the region the assumptions pin (a sign-consistent instantiation), while the
 * residue arithmetic stays symbolic -- so the closed form is correct BY
 * CONSTRUCTION, with no numeric quadrature crosscheck.  A parameter the
 * assumptions leave two-sided unbounded is refused (its sign is undetermined).
 * Without assumptions the classification stays purely numeric (unchanged).
 */
Expr* integrate_residue_try(Expr* f, Expr* x, Expr* a, Expr* b,
                            Expr* assumptions, bool* diverges);

/* `Integrate`ContourResidue[f, {x, a, b}]` builtin.  Strict: NULL on any
 * non-applicable input (no fallback). */
Expr* builtin_integrate_contour_residue(Expr* res);

/* Register the package builtin + attributes + docstring. */
void integrate_residue_init(void);

#endif /* MATHILDA_INTEGRATE_RESIDUE_H */

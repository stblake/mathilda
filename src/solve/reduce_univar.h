/*
 * reduce_univar.h
 *
 * Univariate real sign-diagram engine for `Reduce` (REDUCE_PLAN.md, Phase 2).
 *
 * Solves an arbitrary logical combination (a DNF RForm) of polynomial equations
 * and inequalities in a single variable over the reals.  The method is the
 * classic sign diagram / 1-D cylindrical decomposition:
 *
 *   1. collect the real roots of every polynomial in the formula (via Solve
 *      over Reals) as ordered breakpoints  b_1 < ... < b_m;
 *   2. the reals split into 2m+1 cells -- the open intervals between the
 *      breakpoints and the breakpoints themselves;
 *   3. evaluate the formula's truth at one sample point per cell (a rational
 *      strictly inside an interval, or the breakpoint itself for a point),
 *      using an exact real sign oracle;
 *   4. emit the union of the satisfying cells as a logical combination of
 *      `<`, `<=`, `==`, `!=` atoms (with `Inequality` chains for two-sided
 *      intervals, and `x != r` for isolated excluded points).
 *
 * Reuses `Solve[..., Reals]` for root isolation, `PolynomialQ`/`ReplaceAll`/`N`
 * for the numeric plumbing, and `flint_qqbar_compare` as the exact real-
 * algebraic sign oracle (falling back to exact rational sign for rational
 * samples, so rational-root problems need no FLINT).  Any undecidable sign or
 * ordering makes the engine return NULL, leaving Reduce unevaluated rather than
 * emitting a wrong answer.
 */
#ifndef REDUCE_UNIVAR_H
#define REDUCE_UNIVAR_H

#include "expr.h"
#include "reduce_form.h"

/* Solve the DNF formula `F` in the single real variable `var` over the reals.
 * Returns a freshly-owned Expr (True / False / a logical combination of
 * relations), or NULL to decline (non-polynomial, free parameters, or an
 * undecidable sign/ordering) -- in which case Reduce stays unevaluated. */
Expr* reduce_univar(const RForm* F, const Expr* var, Expr** vars, int nv);

/* Integer-domain variant: the integers satisfying `F`, as an OR of `x == n`
 * atoms, when the solution set is bounded (both unbounded end-cells unsatisfied).
 * Returns True (all integers) / False (none) / an OR, or NULL to decline
 * (unbounded integer set, or an undecidable sign). */
Expr* reduce_univar_integers(const RForm* F, const Expr* var, Expr** vars, int nv);

#endif /* REDUCE_UNIVAR_H */

/*
 * reduce_fm.h
 *
 * Linear real-arithmetic engine for `Reduce` (REDUCE_PLAN.md, Phase 3).
 *
 * Solves a multivariate system of LINEAR equations and inequalities over the
 * reals by Fourier-Motzkin elimination: variables are projected out one at a
 * time (last to first), which both decides feasibility (an eliminated system
 * that reduces to a false constant is unsatisfiable -> False) and yields a
 * triangular description -- bounds on the first variable, then each later
 * variable bounded by expressions in the earlier ones, e.g.
 *
 *   Reduce[x + y < 1 && x > 0 && y > 0, {x, y}, Reals]
 *       ->  0 < x < 1 && 0 < y < 1 - x
 *
 * Equations are handled as paired `<=` constraints and re-detected as `==` in
 * the output.  A disjunction (DNF) is solved conjunct-by-conjunct and OR-ed.
 * Non-linear atoms, disequations (`!=`), rational-function atoms, or non-numeric
 * (parametric) coefficients make the engine decline (return NULL), leaving
 * Reduce for a later phase rather than answering wrongly.
 */
#ifndef REDUCE_FM_H
#define REDUCE_FM_H

#include "expr.h"
#include "reduce_form.h"

/* Solve the DNF formula `F` of linear relations in `vars` (length `nv`) over the
 * reals.  Returns a freshly-owned Expr (True / False / a triangular conjunction,
 * OR-ed over disjuncts), or NULL to decline. */
Expr* reduce_fm(const RForm* F, Expr** vars, int nv);

#endif /* REDUCE_FM_H */

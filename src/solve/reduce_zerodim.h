/*
 * reduce_zerodim.h
 *
 * Zero-dimensional polynomial-system engine, shared by `Reduce` and `Solve`.
 *
 * Where the linear-system engine (reduce_sys) handles affine systems and the
 * CAD engine (reduce_cad) handles positive-dimensional real varieties with
 * rational fibre samples, this engine handles the remaining common case: a
 * conjunction whose polynomial EQUATIONS pin the variety down to a FINITE set
 * of points (a zero-dimensional ideal), optionally carrying polynomial
 * inequality / disequation side relations.
 *
 * Algorithm (general, exact):
 *   1. Split each conjunct into equations E and side relations O (<, <=, !=).
 *   2. Solve E over the complexes with the existing polynomial-system solver
 *      (reused via an internal Solve[] call).  A finite, fully-determined
 *      solution set means E is zero-dimensional; a parametric / underdetermined
 *      answer means it is not, and the engine DECLINES.
 *   3. Keep each solution branch only if every side relation holds there and --
 *      over the Reals -- every coordinate is real, decided EXACTLY with the
 *      FLINT qqbar algebraic-number oracle (sign / equality / realness).  An
 *      undecidable test makes the engine DECLINE (sound over complete).
 *
 * This is complete for zero-dimensional systems (the finite solution set is
 * enumerated and filtered exactly) and sound (any gap yields a decline, leaving
 * the input unevaluated).
 */
#ifndef REDUCE_ZERODIM_H
#define REDUCE_ZERODIM_H

#include "expr.h"
#include "reduce_form.h"
#include "reduce_opts.h"
#include "solvepoly.h"   /* SolvePolyOpts */
#include <stdbool.h>

/* Reduce engine.  `F` is the DNF form, `vars`/`nv` the ambient variables,
 * `reals` selects the Reals (true) or Complexes (false) domain.  Returns the
 * And/Or solution Expr (freshly owned), or NULL to decline. */
Expr* reduce_zerodim(const RForm* F, Expr** vars, int nv, bool reals,
                     const ReduceOpts* opts);

/* Solve engine.  `expr` is the evaluated first argument (an And of relations,
 * or a bare relation), `vars` a symbol or List of symbols, `dom` the domain
 * symbol or NULL (Complexes by default, but an ordering inequality forces the
 * Reals).  Returns a List of solution rule-lists (List[List[Rule[var,val]]]) --
 * the empty List for a provably empty set -- or NULL to decline.  Declines when
 * the conjunction carries no inequality / disequation side relation (leaving a
 * pure equation system to Solve's own specialists). */
Expr* reduce_zerodim_solve(const Expr* expr, const Expr* vars, const Expr* dom,
                           const SolvePolyOpts* poly);

#endif /* REDUCE_ZERODIM_H */

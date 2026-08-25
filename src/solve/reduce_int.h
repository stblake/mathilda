/*
 * reduce_int.h
 *
 * Integer / Rationals domain for `Reduce` (REDUCE_PLAN.md, Phase 5).
 *
 * A thin adapter over the existing `Solve[..., Integers|Rationals]` engine: its
 * solution list (finite tuples, or parametric families with C[k] generators) is
 * reformatted into Reduce's logical form -- an OR of conjunctions of `var ==
 * value` atoms, with `Element[C[k], dom]` added for each generated parameter.
 *
 * When Solve declines (e.g. a pure bounded inequality system, which the integer
 * engine only handles alongside an equation), the univariate case falls back to
 * bounded integer enumeration over the real sign diagram
 * (reduce_univar_integers).
 */
#ifndef REDUCE_INT_H
#define REDUCE_INT_H

#include "expr.h"
#include "reduce_form.h"
#include "reduce_opts.h"

/* Solve `expr` (the original input relation) for `vars_expr` over `dom`
 * (Integers or Rationals) and return the logical form, or NULL to decline.
 * `F`/`vars`/`nv` back the univariate enumeration fallback.  `opts` supplies
 * GeneratedParameters (the head used for free parameters). */
Expr* reduce_integers(const Expr* expr, const Expr* vars_expr, const RForm* F,
                      Expr** vars, int nv, const Expr* dom,
                      const ReduceOpts* opts);

/* Modulus -> p residue enumeration.  Solves `expr` for `vars_expr` over Z/pZ
 * (p = opts->modulus) via Solve's modular engine and reformats the finite
 * solution list into an Or of `var == r` equations.  Returns NULL (Reduce
 * stays unevaluated) when the modulus is symbolic / out of range or the
 * statement is not solvable modularly. */
Expr* reduce_modular(const Expr* expr, const Expr* vars_expr,
                     const ReduceOpts* opts);

#endif /* REDUCE_INT_H */

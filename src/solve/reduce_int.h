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

/* Solve `expr` (the original input relation) for `vars_expr` over `dom`
 * (Integers or Rationals) and return the logical form, or NULL to decline.
 * `F`/`vars`/`nv` back the univariate enumeration fallback. */
Expr* reduce_integers(const Expr* expr, const Expr* vars_expr, const RForm* F,
                      Expr** vars, int nv, const Expr* dom);

#endif /* REDUCE_INT_H */

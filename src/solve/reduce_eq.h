/*
 * reduce_eq.h
 *
 * The complete univariate equation engine for `Reduce` (docs/notes/REDUCE_PLAN.md, Phase 1).
 *
 * `Solve` returns only the GENERIC solution of a polynomial equation; `Reduce`
 * must report the COMPLETE solution set, including every case where a leading
 * coefficient vanishes.  reduce_eq_univariate builds that as a DNF tree:
 *
 *   solve(p in x):
 *     lc = leading coeff of p in x
 *     lc a nonzero constant  ->  (generic solutions of p)
 *     otherwise              ->  (lc != 0 && generic solutions of p)
 *                             || (lc == 0 && solve(p with leading term dropped))
 *
 * The generic solutions come straight from `Solve` (reused); the vanishing-
 * coefficient recursion is the piece Solve does not provide.  This mirrors the
 * coefficient-vanishing recursion of `SolveAlways` (src/solve/solvealways.c).
 */
#ifndef REDUCE_EQ_H
#define REDUCE_EQ_H

#include "expr.h"
#include "reduce_form.h"
#include "reduce_opts.h"

/* Complete solution set of a single univariate polynomial equation `poly == 0`
 * in `var` over the complex numbers, as a DNF RForm.  `vars`/`nv` are the
 * ambient reduce variables (for atom classification).  `opts` forwards the
 * Cubics / Quartics radical flags onto the internal Solve calls.  Returns a
 * freshly-owned RForm; sets *ok = false (leaving the form meaningless, to be
 * freed) when the generic solver declines or hands back an unexpected shape,
 * so the caller can leave the whole Reduce unevaluated. */
RForm* reduce_eq_univariate(const Expr* poly, const Expr* var,
                            Expr** vars, int nv, bool* ok,
                            const ReduceOpts* opts);

/* Complete solution set of a single univariate TRANSCENDENTAL equation
 * `poly == 0` in `var` (Log/Exp/inverse-trig over `var`), rendered directly as
 * a logical-formula Expr* by re-entering `Solve` (which now combines multi-log
 * residuals and inverts Exp kernels). The rule-list `{{var -> v}, ...}` becomes
 * an Or of  `var == v`  disjuncts; a `ConditionalExpression[v, cond]` solution
 * (periodic `Element[C,Integers]` or the log/inverse-trig strip
 * `-Pi < Im[C] <= Pi`) becomes `cond && var == v`.
 *
 * Unlike reduce_eq_univariate this returns a raw Expr* rather than an RForm,
 * because those branch conditions are not polynomial atoms. Returns NULL when
 * Solve declines or hands back an unexpected shape, so the caller falls back to
 * the polynomial / echo path. An empty solution set yields False; a `{{}}`
 * (all-values) row yields True. */
Expr* reduce_eq_transcendental(const Expr* poly, const Expr* var,
                               const Expr* dom, const ReduceOpts* opts);

/* Realness of a solution value `v` over Reals, with every free parameter real:
 * 1 real, 0 never real, 2 real exactly when *cond_out (owned, evaluated; a
 * condition on the parameters such as -1 <= a <= 1 for ArcSin[a]) holds, -1
 * undecidable.  *cond_out is set only for the result 2. */
int reduce_value_reality(const Expr* v, Expr** cond_out);

/* Over Reals, reduce a Solve rule list `{{var -> v}, ...}` (consumed) to its
 * real members: a periodic family a + b C[k] is kept whole when its period and
 * offset are real, collapsed to its C[k] = 0 member when the period is not
 * real, and dropped when it has no real member.  A value that is real only for
 * some values of the (real) parameters -- ArcSin[a], Log[a] -- is kept with
 * its exact reality condition conjoined onto the row's ConditionalExpression.
 * Returns NULL (the caller declines) when the reality of a value can be
 * neither decided nor expressed as a condition. */
Expr* reduce_collapse_reals_families(Expr* sols);

#endif /* REDUCE_EQ_H */

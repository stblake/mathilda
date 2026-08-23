/*
 * reduce_eq.h
 *
 * The complete univariate equation engine for `Reduce` (REDUCE_PLAN.md, Phase 1).
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

/* Complete solution set of a single univariate polynomial equation `poly == 0`
 * in `var` over the complex numbers, as a DNF RForm.  `vars`/`nv` are the
 * ambient reduce variables (for atom classification).  Returns a freshly-owned
 * RForm; sets *ok = false (leaving the form meaningless, to be freed) when the
 * generic solver declines or hands back an unexpected shape, so the caller can
 * leave the whole Reduce unevaluated. */
RForm* reduce_eq_univariate(const Expr* poly, const Expr* var,
                            Expr** vars, int nv, bool* ok);

#endif /* REDUCE_EQ_H */

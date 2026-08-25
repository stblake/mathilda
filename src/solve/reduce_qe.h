/*
 * reduce_qe.h
 *
 * Quantifier elimination for `Reduce` (REDUCE_PLAN.md, Phase 7): the `Exists`,
 * `ForAll` and `Resolve` heads.  A quantified real statement is turned into an
 * equivalent quantifier-free formula over the remaining free variables (or
 * True / False when none remain), by Cylindrical Algebraic Decomposition with
 * the quantified variables projected out first.
 *
 * Three cases by the number of free variables (nfree):
 *   - nfree == 0 (fully quantified): a real-closed-field DECISION procedure --
 *     Exists is `Reduce[phi, bound, Reals] =!= False`, ForAll is `... === True`
 *     -- reusing the whole existing engine.
 *   - nfree == 1 (parametric, single quantifier block): reduce_cad_qe, the CAD
 *     seam that folds each free-variable cell's bound subtree to an
 *     Exists/ForAll verdict and emits the 1-D sign diagram.
 *   - nfree >= 2, alternating quantifiers, or an explicit non-Reals domain:
 *     decline (NULL) in v1 (needs nested QE emission / Phase 6b fibres).
 *
 * Soundness invariant (shared with the rest of Reduce): an undecidable sign, a
 * non-rational fibre, an unsupported construct, or an out-of-scope case makes the
 * engine return NULL, leaving the input unevaluated rather than guessed.
 */
#ifndef REDUCE_QE_H
#define REDUCE_QE_H

#include "expr.h"

/* Register Exists / ForAll / Resolve.  Called from reduce_init. */
void  reduce_qe_init(void);

/* `Resolve[qexpr]` / `Resolve[qexpr, dom]` -- eliminate the quantifiers in
 * qexpr.  Returns NULL (unevaluated) for a non-quantified qexpr. */
Expr* builtin_resolve(Expr* res);

/* Inert quantifier wrappers: Exists / ForAll carry a binding for Reduce and
 * Resolve to eliminate and do not evaluate on their own (return NULL). */
Expr* builtin_exists(Expr* res);
Expr* builtin_forall(Expr* res);

/* Eliminate the quantifiers in `qexpr` (a top-level Exists/ForAll) over `dom`
 * (borrowed; NULL selects the default Reals; anything but Reals declines in v1).
 * Returns a freshly-owned quantifier-free Expr, or NULL to decline.  Used by
 * Resolve and by Reduce's Exists/ForAll front-end peel. */
Expr* reduce_qe_dispatch(const Expr* qexpr, const Expr* dom);

#endif /* REDUCE_QE_H */

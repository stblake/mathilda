#ifndef MATHILDA_SIMP_LOG_H
#define MATHILDA_SIMP_LOG_H

#include "expr.h"
#include "simp.h"

/*
 * simp_log -- logarithm-specific simplifications for Simplify.
 *
 * Implements two orthogonal primitives that together close the common
 * Log identities Simplify would otherwise miss:
 *
 *   Pass A -- Prime decomposition. Every Log[r] subexpression whose
 *     argument is a positive rational r = p/q is replaced by the integer
 *     linear combination Sum e_i Log[p_i] - Sum f_j Log[q_j] obtained by
 *     prime-factoring p and q. Also recognises Log[Power[r, k]] for
 *     positive-rational r, rewriting as k * (decomposed Log[r]). Sound
 *     unconditionally; closes Log[4] -> 2 Log[2], Log[72] -> 3 Log[2] +
 *     2 Log[3], Log[2/3] -> Log[2] - Log[3], Log[Sqrt[12]] ->
 *     Log[2] + (1/2) Log[3], and similar.
 *
 *   Pass B -- Linear-combination-of-logs fuser. For every Plus
 *     subexpression containing >=2 terms of the form c_i * Log[a_i],
 *     fuse the log block into a single Log[ Product a_i ^ c_i ], let
 *     the evaluator + Together cancel inside the product, and take the
 *     fused form only when it is strictly simpler. Gated on positivity
 *     of every a_i (consults the AssumeCtx) when the product is
 *     symbolic; constant-collapse fusions (result is a leaf or near-
 *     leaf) are taken unconditionally as they cannot introduce new
 *     branch-cut hazards.
 *
 * Both passes are run to a fixed point (bounded). Pass A is idempotent
 * after a single application but Pass B can expose new
 * decomposable-rational arguments for Pass A to eat (e.g. the fused
 * product evaluates to a positive rational).
 *
 * Returns NULL if the input is unchanged. Otherwise returns a freshly
 * owned, evaluated expression tree. `ctx` may be NULL (treated as an
 * empty assumption context).
 */
Expr* simp_log_apply(const Expr* e, const AssumeCtx* ctx);

/*
 * simp_log_fuse_all -- unconditional linear-combination-of-logs fuser, for
 * equation solving (not Simplify).
 *
 * Given a top-level Plus e = Sum_i c_i Log[a_i] + rest, fuse the Log terms
 * whose argument contains `var` (or ALL Log terms when `var` is NULL) into a
 * single Log[ Product a_i ^ c_i ], applying Together+Cancel to the fused
 * argument so algebraic cancellations surface (e.g.
 * Log[t^2-x^2] - Log[t-x] -> Log[t+x]). Returns Plus[fused_log, rest...],
 * freshly owned; returns NULL when fewer than two matching Log terms are
 * present (nothing to fuse).
 *
 * Unlike Pass B (try_fuse_plus / simp_log_apply), this takes the fusion
 * UNCONDITIONALLY -- no positivity gate and no leaf-count improvement test --
 * because collapsing a multi-log equation to one invertible Log head is what
 * lets the inverse-function solver peel it, even when the fused form is not
 * smaller by leaf count. The fusion is principal-branch (as is Log[a]+Log[b]
 * -> Log[a b]); callers that surface it to the user emit Solve::ifun.
 *
 * `var`, when non-NULL, must be an EXPR_SYMBOL.
 */
Expr* simp_log_fuse_all(const Expr* e, const Expr* var);

#endif /* MATHILDA_SIMP_LOG_H */

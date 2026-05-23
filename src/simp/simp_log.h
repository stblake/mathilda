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

#endif /* MATHILDA_SIMP_LOG_H */

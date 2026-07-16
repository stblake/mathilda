/* cherry_dilog.h — dilogarithm integration (Cherry degree-2 Sigma-decomposition).
 *
 * R(x) Log[w] -> Log-Log products + PolyLog[2, g] (dilogarithms), matched in the
 * log tower with dilog arguments g the degree-1 Sigma-decomposition candidates
 * between the rational roots of the linear factors.  Generalises rt_try_dilog to
 * the Log-Log + PolyLog[2] answer form (INT Log[x]/(1+x) = Log[x]Log[1+x] +
 * PolyLog[2,-x]).  Returns a fresh, PowerExpand-diff-back-verified antiderivative
 * or NULL.  Defined in cherry_dilog.c.
 */

#ifndef MATHILDA_CHERRY_DILOG_H
#define MATHILDA_CHERRY_DILOG_H

#include "expr.h"

Expr* rt_cherry_dilog(Expr* f, Expr* x);

#endif /* MATHILDA_CHERRY_DILOG_H */

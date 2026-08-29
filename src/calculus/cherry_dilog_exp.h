/* cherry_dilog_exp.h — exponential-tower dilogarithm integration (Cherry).
 *
 * The exponential-tower mirror of rt_cherry_dilog (cherry_dilog.c): integrands
 * that are a rational-in-theta combination of the weight-1 logs {x} u {Log[theta
 * - rho_k]}, with theta = E^(c x), whose antiderivative is a combination of
 * Log-Log products, x-weighted logs, and dilogarithms PolyLog[2, g] with g a
 * Moebius function of theta.  Here x = Log[theta]/c is the "root at 0" tower-log
 * and Log[1 + E^x] = Log[theta - (-1)] is the tower-log for rho = -1, so both the
 * rational-times-x forms (x/(E^x - 1)) and the outer-log forms (Log[1 + E^x])
 * are the SAME input form, matched natively in the tower — no substitution.
 *
 *   INT x/(E^x - 1) dx = x Log[1 - E^-x] - PolyLog[2, E^-x]
 *   INT Log[1 + E^x] dx = -PolyLog[2, -E^x]
 *
 * Returns a fresh, PowerExpand-diff-back-verified antiderivative or NULL.
 * Defined in cherry_dilog_exp.c.
 */

#ifndef MATHILDA_CHERRY_DILOG_EXP_H
#define MATHILDA_CHERRY_DILOG_EXP_H

#include "expr.h"

Expr* rt_cherry_dilog_exp(Expr* f, Expr* x);

/* Integrate`Cherry`DilogExp[f, x] — direct debuggable surface for the engine
 * (bypasses the cascade), used by tests and benchmarks.  Registered by
 * cherry_builtins_init (cherry_driver.c). */
Expr* builtin_cherry_dilog_exp(Expr* res);

#endif /* MATHILDA_CHERRY_DILOG_EXP_H */

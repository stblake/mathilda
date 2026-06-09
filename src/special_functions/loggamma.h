/* Mathilda -- LogGamma[z], the log-gamma function log(Gamma(z)).
 *
 * LogGamma is the analytic continuation of log(Gamma(z)) with a single branch
 * cut along the negative real axis -- it is NOT identical to Log[Gamma[z]],
 * which inherits the (more intricate) branch structure of Log. Unlike Gamma,
 * LogGamma stays finite and numerically well-behaved for very large arguments
 * (where Gamma itself overflows), so it is the right primitive for ratios of
 * factorials and combinatorial asymptotics.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#ifndef MATHILDA_LOGGAMMA_H
#define MATHILDA_LOGGAMMA_H

#include "expr.h"

/* Builtin entry point: LogGamma[z]. Takes ownership of res per the builtin
 * contract (returns a new Expr* on success, NULL to stay unevaluated). */
Expr* builtin_loggamma(Expr* res);

/* Register the builtin, its attributes and docstring. */
void loggamma_init(void);

#endif /* MATHILDA_LOGGAMMA_H */

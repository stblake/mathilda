#ifndef EXPINTEGRALE_H
#define EXPINTEGRALE_H

#include "expr.h"

/* ExpIntegralE[n, z] -- the generalized exponential integral
 *   E_n(z) = Integrate[E^(-z t) / t^n, {t, 1, Infinity}],   Re z > 0.
 * E_0(z) = e^(-z)/z; E_n'(z) = -E_(n-1)(z); E_1(z) = Gamma[0, z]. */
Expr* builtin_expintegrale(Expr* res);
void  expintegrale_init(void);

/* Machine E_n(x) for a nonnegative integer order n and real x > 0 (the classic
 * continued-fraction / series split).  Returns false (defer) outside that
 * domain.  Shared with the NDArray kernel. */
bool expintegrale_machine(long n, double x, double* out);

#endif /* EXPINTEGRALE_H */

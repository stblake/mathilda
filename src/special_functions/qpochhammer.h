/* Mathilda -- QPochhammer, the q-Pochhammer symbol (q-shifted factorial).
 *
 *   QPochhammer[a, q, n] = prod_{k=0}^{n-1} (1 - a q^k)        (finite)
 *   QPochhammer[a, q]     = prod_{k=0}^{Infinity} (1 - a q^k)  ((a;q)_inf, |q|<1)
 *
 * The finite form is exact/symbolic for a non-negative integer n (it expands
 * to a product the evaluator reduces, so machine and MPFR N come for free);
 * the infinite form evaluates numerically for machine-real a, q with |q| < 1.
 */
#ifndef QPOCHHAMMER_H
#define QPOCHHAMMER_H

#include "expr.h"

Expr* builtin_qpochhammer(Expr* res);
void qpochhammer_init(void);

#endif /* QPOCHHAMMER_H */

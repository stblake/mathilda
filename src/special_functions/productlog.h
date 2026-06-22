#ifndef PRODUCTLOG_H
#define PRODUCTLOG_H

#include "expr.h"

/* ProductLog -- the Lambert W function.
 *
 *   ProductLog[z]     principal branch W_0(z): the solution w of z == w e^w
 *                     whose value is real for z >= -1/e.
 *   ProductLog[k, z]  the k-th branch W_k(z) (k any integer, k == 0 principal).
 *
 * Numeric evaluation (real/complex, machine and arbitrary MPFR precision) runs
 * through a single complex-MPFR Halley core; exact special values short-circuit
 * it. Attributes: Listable, NumericFunction, Protected, ReadProtected.
 * D[ProductLog[z], z] = ProductLog[z] / (z (1 + ProductLog[z])). */
Expr* builtin_productlog(Expr* res);
void  productlog_init(void);

#endif /* PRODUCTLOG_H */

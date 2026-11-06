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
 * it. Attributes: Listable, NumericFunction, Protected.
 * D[ProductLog[z], z] = ProductLog[z] / (z (1 + ProductLog[z])). */
Expr* builtin_productlog(Expr* res);
void  productlog_init(void);
#ifdef USE_MPFR
#include <mpfr.h>
/* Real principal-branch W at rop's precision (mpfr unary signature); NaN when
 * x is outside (-1/e + 2^-10, oo) or not a number. */
int productlog_mpfr_real(mpfr_ptr rop, mpfr_srcptr x, mpfr_rnd_t rnd);
#endif

#endif /* PRODUCTLOG_H */

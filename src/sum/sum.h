#ifndef SUM_H
#define SUM_H

#include "expr.h"

/*
 * Sum[f, {i, imax}]                evaluate sum_{i=1}^{imax} f
 * Sum[f, {i, imin, imax}]          sum_{i=imin}^{imax} f
 * Sum[f, {i, imin, imax, di}]      step di
 * Sum[f, {i, {i1, i2, ...}}]       successive values
 * Sum[f, {i,...}, {j,...}, ...]    multiple (nested) sums
 * Sum[f, i]                        indefinite sum (antidifference)
 *
 * Implemented as a dispatcher (this file) that handles iteration and a
 * Method cascade over the context-qualified sub-algorithms
 * Sum`Polynomial, Sum`Geometric, Sum`Gosper, ... -- mirroring how Integrate
 * dispatches to Integrate`BronsteinRational etc.
 */
Expr* builtin_sum(Expr* res);

/*
 * Convergence of an infinite series / product (sum_converge.c).  f is the
 * evaluated term at the symbolic index k (the caller shadows k); only the
 * tail k -> Infinity is examined.  UNKNOWN never means "diverges".
 */
typedef enum {
    SUM_CONV_UNKNOWN = 0,
    SUM_CONV_YES,           /* the series converges                  */
    SUM_CONV_NO             /* it diverges (Sum::div)                */
} SumConvergence;

typedef enum {
    PROD_CONV_UNKNOWN = 0,
    PROD_CONV_YES,          /* converges to a nonzero limit          */
    PROD_CONV_DIV,          /* tends to Infinity or has no limit     */
    PROD_CONV_ZERO          /* tends to 0 (the product IS 0)         */
} ProductConvergence;

SumConvergence sum_series_convergence(Expr* f, Expr* k);
ProductConvergence product_convergence(Expr* f, Expr* k);

void sum_init(void);

#endif

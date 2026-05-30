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

void sum_init(void);

#endif

#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "expr.h"

/* LegendreP -- Legendre polynomials P_n(x) and the associated Legendre
 * functions of the first kind P_n^m(x).
 *
 *   LegendreP[n, x]        Legendre polynomial / function P_n(x).
 *   LegendreP[n, m, x]     associated Legendre function P_n^m(x) (type 1).
 *   LegendreP[n, m, a, x]  Legendre function of type a (a in {1, 2, 3};
 *                          a == 1 is the default).
 *
 * Exact integer n give the explicit polynomial in x (exact rational
 * coefficients). Non-integer order with an inexact argument is evaluated
 * numerically through the Gauss series 2F1(-n, n+1; 1; (1-x)/2) at machine or
 * arbitrary MPFR precision (real and complex). The associated functions for
 * integer n and m >= 0 use the Rodrigues derivative form for type 1 and the
 * regularized 2F1 form for types 2 and 3.
 *
 * Attributes: Listable, NumericFunction, Protected. */
Expr* builtin_legendre_p(Expr* res);
void  legendre_init(void);

#endif /* LEGENDRE_H */

#ifndef SUM_EULER_INTERNAL_H
#define SUM_EULER_INTERNAL_H

#include "expr.h"

/*
 * Shared builders / depth-<=2 reducers used by both sum_euler.c (linear Euler
 * sums) and sum_euler_nonlinear.c (nonlinear Euler sums via MZV reduction).
 * All return owned, still-unevaluated Expr* trees in Riemann zeta values.
 */

/* Zeta[n] node. */
Expr* eu_zeta(int n);
/* 1/2 as Power[2,-1]. */
Expr* eu_half(void);
/* Times[a,b] / Plus[a,b] (adopt a,b). */
Expr* eu_times2(Expr* a, Expr* b);
Expr* eu_plus2(Expr* a, Expr* b);
/* Exact integer binomial C(n,k). */
int64_t eu_binom(int n, int k);

/* Euler's Sum_{k>=1} H_k / k^q  (q >= 2), unevaluated. */
Expr* euler_order1(int q);
/* Double zeta Z(s,t) = zeta(s,t) for ODD s (Borwein-Borwein-Girgensohn). */
Expr* euler_Z_odd(int s, int t);

#endif /* SUM_EULER_INTERNAL_H */

/* Mathilda -- JacobiP, the Jacobi polynomial P_n^(a,b)(x).
 *
 *   JacobiP[n, a, b, x]   gives the Jacobi polynomial P_n^(a,b)(x), orthogonal
 *                         on [-1, 1] with weight (1-x)^a (1+x)^b.
 *
 * See jacobi_p.c for the evaluation strategy (exact integer-degree polynomial,
 * stable numeric recurrence at machine and arbitrary precision, and the Gauss
 * hypergeometric form for a non-integer degree).
 */
#ifndef MATHILDA_JACOBI_P_H
#define MATHILDA_JACOBI_P_H

#include "expr.h"

/* Builtin entry point. Takes ownership of `res` per the builtin contract. */
Expr* builtin_jacobi_p(Expr* res);

/* Register JacobiP, its attributes and docstring. Called from core_init(). */
void jacobi_p_init(void);

#endif /* MATHILDA_JACOBI_P_H */

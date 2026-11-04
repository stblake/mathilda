/*
 * matfunc.h -- the matrix functions: MatrixExp, MatrixLog, MatrixFunction and
 * the non-integer exponents of MatrixPower.
 *
 * Standard builtin ownership contract (SPEC.md section 4): each takes `res`,
 * returns a fresh Expr* on success or NULL to leave the call unevaluated, and
 * never frees `res`.
 */
#ifndef MATHILDA_LINALG_MATFUNC_H
#define MATHILDA_LINALG_MATFUNC_H

#include "expr.h"

Expr* builtin_matrixexp(Expr* res);        /* MatrixExp[m], MatrixExp[m, v]   */
Expr* builtin_matrixlog(Expr* res);        /* MatrixLog[m]                    */
Expr* builtin_matrixfunction(Expr* res);   /* MatrixFunction[f, m]            */

/* MatrixPower[m, s] / MatrixPower[m, s, v] for an exponent s that is not an
 * integer (a rational, real, complex or symbolic s): the principal power,
 * through the same machinery as MatrixFunction.  Called by builtin_matrixpower
 * (matpow.c) after it has validated m; the integer path never comes here. */
Expr* matfunc_matrixpower(Expr* res);

void matfunc_init(void);

#endif /* MATHILDA_LINALG_MATFUNC_H */

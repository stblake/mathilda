#ifndef LINALG_H
#define LINALG_H

#include <stddef.h>
#include <stdint.h>
#include "expr.h"

Expr* builtin_dot(Expr* res);
Expr* builtin_det(Expr* res);
Expr* builtin_cross(Expr* res);
Expr* builtin_norm(Expr* res);
Expr* builtin_tr(Expr* res);
Expr* builtin_identitymatrix(Expr* res);
Expr* builtin_diagonalmatrix(Expr* res);
Expr* builtin_inverse(Expr* res);
Expr* builtin_matrixpower(Expr* res);
Expr* builtin_eigenvalues(Expr* res);
Expr* builtin_eigenvectors(Expr* res);
void linalg_init(void);

/* Helpers exposed for use by matsol.c (RowReduce / LinearSolve). */
int   get_tensor_dims(Expr* e, int64_t* dims);
void  flatten_tensor(Expr* e, Expr** flat, size_t* idx);
Expr* laplace_det(Expr** flat, int original_n, int n, int row, int* cols);
Expr* exact_div_wrapper(Expr* num, Expr* den);

#endif // LINALG_H
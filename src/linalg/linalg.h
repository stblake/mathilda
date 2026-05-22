#ifndef LINALG_H
#define LINALG_H

#include <stddef.h>
#include <stdint.h>
#include "expr.h"

/* Built-in entry points registered by linalg_init().  Per-builtin
 * implementations live in their own translation units inside
 * src/linalg/ (dot.c, det.c, cross.c, ...). */
Expr* builtin_dot(Expr* res);
Expr* builtin_det(Expr* res);
Expr* builtin_cross(Expr* res);
Expr* builtin_norm(Expr* res);
Expr* builtin_normalize(Expr* res);
Expr* builtin_tr(Expr* res);
Expr* builtin_identitymatrix(Expr* res);
Expr* builtin_diagonalmatrix(Expr* res);
/* Inverse and PseudoInverse live in src/matinv.h. */
Expr* builtin_matrixpower(Expr* res);
/* Eigenvalues / Eigenvectors live in src/mateigen.h. */
Expr* builtin_positive_definite_matrix_q(Expr* res);
Expr* builtin_negative_definite_matrix_q(Expr* res);
void linalg_init(void);

/* Helpers exposed for use by matsol.c (RowReduce / LinearSolve) and
 * mateigen.c (Eigenvalues / Eigenvectors). */
int   get_tensor_dims(Expr* e, int64_t* dims);
void  flatten_tensor(Expr* e, Expr** flat, size_t* idx);
Expr* laplace_det(Expr** flat, int original_n, int n, int row, int* cols);
Expr* exact_div_wrapper(Expr* num, Expr* den);
/* dot2: internal tensor-contraction helper used by Dot, MatrixPower, and the
 * eigenvalue solver.  Returns a freshly-allocated Expr*, or NULL on shape
 * mismatch (in which case *error_printed is set when a diagnostic was
 * already emitted).  Callers own the result. */
#ifndef __cplusplus
#include <stdbool.h>
#endif
Expr* dot2(Expr* a, Expr* b, bool* error_printed);

#endif // LINALG_H

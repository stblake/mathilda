#ifndef MATRIXQ_H
#define MATRIXQ_H

#include "expr.h"

Expr* builtin_hermitian_matrix_q(Expr* res);
Expr* builtin_symmetric_matrix_q(Expr* res);
Expr* builtin_square_matrix_q(Expr* res);
Expr* builtin_diagonal_matrix_q(Expr* res);
Expr* builtin_upper_triangular_matrix_q(Expr* res);

#endif /* MATRIXQ_H */

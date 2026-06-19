#ifndef LIST_H
#define LIST_H

#include "expr.h"

Expr* builtin_table(Expr* res);
Expr* builtin_range(Expr* res);
Expr* builtin_array(Expr* res);
Expr* builtin_take(Expr* res);
Expr* builtin_drop(Expr* res);
Expr* builtin_flatten(Expr* res);
Expr* builtin_partition(Expr* res);
Expr* builtin_rotateleft(Expr* res);
Expr* builtin_rotateright(Expr* res);
Expr* builtin_reverse(Expr* res);
Expr* builtin_rescale(Expr* res);
Expr* builtin_transpose(Expr* res);
Expr* builtin_conjugate_transpose(Expr* res);
Expr* builtin_tally(Expr* res);
Expr* builtin_union(Expr* res);
Expr* builtin_deleteduplicates(Expr* res);
Expr* builtin_split(Expr* res);
Expr* builtin_total(Expr* res);
Expr* builtin_accumulate(Expr* res);
Expr* builtin_differences(Expr* res);
Expr* builtin_ratios(Expr* res);
Expr* builtin_commonest(Expr* res);
Expr* builtin_min(Expr* res);
Expr* builtin_max(Expr* res);
Expr* builtin_listq(Expr* res);
Expr* builtin_vectorq(Expr* res);
Expr* builtin_matrixq(Expr* res);
Expr* builtin_hermitian_matrix_q(Expr* res);
Expr* builtin_symmetric_matrix_q(Expr* res);
Expr* builtin_square_matrix_q(Expr* res);
Expr* builtin_diagonal_matrix_q(Expr* res);
Expr* builtin_upper_triangular_matrix_q(Expr* res);

void list_init(void);

Expr* builtin_join(Expr* res);

#endif // LIST_H
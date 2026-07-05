#ifndef SORT_H
#define SORT_H

#include "expr.h"

Expr* builtin_sort(Expr* res);
Expr* builtin_sort_by(Expr* res);
Expr* builtin_maximal_by(Expr* res);
Expr* builtin_minimal_by(Expr* res);
Expr* builtin_take_largest(Expr* res);
Expr* builtin_take_smallest(Expr* res);
Expr* builtin_take_largest_by(Expr* res);
Expr* builtin_take_smallest_by(Expr* res);
Expr* builtin_reverse_sort(Expr* res);
Expr* builtin_reverse_sort_by(Expr* res);
Expr* builtin_orderedq(Expr* res);

#endif // SORT_H

#ifndef PATTERNS_H
#define PATTERNS_H

#include "expr.h"

Expr* builtin_cases(Expr* res);
Expr* builtin_first_case(Expr* res);
Expr* builtin_first_position(Expr* res);
Expr* builtin_delete_missing(Expr* res);
Expr* builtin_delete_cases(Expr* res);
Expr* builtin_position(Expr* res);
Expr* builtin_count(Expr* res);
Expr* builtin_memberq(Expr* res);

void patterns_init(void);

/* LinkedList walk for Count / MemberQ over the elements at `level` (typed-slot
 * fast paths); the LinkedList ds_ops row's descend slot. */
int64_t ll_descend_expr(Expr* x, int64_t level, int64_t min_l, int64_t max_l,
                        bool heads, Expr* pattern, bool first_only);

#endif

#ifndef PATTERNS_H
#define PATTERNS_H

#include "expr.h"

Expr* builtin_cases(Expr* res);
Expr* builtin_first_case(Expr* res);
Expr* builtin_delete_missing(Expr* res);
Expr* builtin_delete_cases(Expr* res);
Expr* builtin_position(Expr* res);
Expr* builtin_count(Expr* res);
Expr* builtin_memberq(Expr* res);

void patterns_init(void);

#endif

#ifndef DATETIME_H
#define DATETIME_H

#include "expr.h"

Expr* builtin_timing(Expr* res);
Expr* builtin_repeated_timing(Expr* res);
Expr* builtin_absolute_time(Expr* res);

void datetime_init(void);

#endif // DATETIME_H
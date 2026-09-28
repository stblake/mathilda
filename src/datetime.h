#ifndef DATETIME_H
#define DATETIME_H

#include "expr.h"

Expr* builtin_timing(Expr* res);
Expr* builtin_absolute_timing(Expr* res);
Expr* builtin_repeated_timing(Expr* res);
Expr* builtin_absolute_time(Expr* res);
Expr* builtin_pause(Expr* res);
Expr* builtin_session_time(Expr* res);
Expr* builtin_time_used(Expr* res);

void datetime_init(void);

#endif // DATETIME_H
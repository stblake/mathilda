#ifndef LOGEXP_H
#define LOGEXP_H

#include "expr.h"

Expr* builtin_log(Expr* res);
Expr* builtin_exp(Expr* res);
Expr* builtin_log10(Expr* res);
Expr* builtin_log2(Expr* res);

void logexp_init(void);

#endif // LOGEXP_H

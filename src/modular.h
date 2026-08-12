#ifndef MODULAR_H
#define MODULAR_H

#include "expr.h"

Expr* builtin_module(Expr* res);
Expr* builtin_block(Expr* res);
Expr* builtin_with(Expr* res);
Expr* builtin_unique(Expr* res);

void modular_init(void);

#endif // MODULAR_H

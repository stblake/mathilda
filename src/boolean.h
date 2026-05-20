#ifndef BOOLEAN_H
#define BOOLEAN_H

#include "expr.h"

Expr* builtin_not(Expr* res);
Expr* builtin_and(Expr* res);
Expr* builtin_or(Expr* res);
Expr* builtin_boole(Expr* res);

// Initialize boolean builtins
void boolean_init(void);

#endif // BOOLEAN_H

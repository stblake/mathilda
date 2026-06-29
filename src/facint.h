#ifndef FACINT_H
#define FACINT_H

#include "expr.h"

Expr* builtin_primeq(Expr* res);
Expr* builtin_nextprime(Expr* res);
Expr* builtin_factorinteger(Expr* res);
Expr* builtin_eulerphi(Expr* res);

void facint_init(void);

#endif

#ifndef BESSEL_H
#define BESSEL_H

#include "expr.h"

Expr* builtin_besselj(Expr* res);
Expr* builtin_besselk(Expr* res);
Expr* builtin_besseli(Expr* res);
Expr* builtin_bessely(Expr* res);
void  bessel_init(void);

#endif /* BESSEL_H */

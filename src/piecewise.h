#ifndef PIECEWISE_H
#define PIECEWISE_H

#include "expr.h"

Expr* builtin_floor(Expr* res);
Expr* builtin_ceiling(Expr* res);
Expr* builtin_round(Expr* res);
Expr* builtin_integerpart(Expr* res);
Expr* builtin_fractionalpart(Expr* res);
Expr* builtin_unitstep(Expr* res);

void piecewise_init(void);

#endif // PIECEWISE_H

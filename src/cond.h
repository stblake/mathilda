#ifndef COND_H
#define COND_H

#include "expr.h"

Expr* builtin_if(Expr* res);
Expr* builtin_which(Expr* res);
Expr* builtin_switch(Expr* res);
Expr* builtin_trueq(Expr* res);
Expr* builtin_piecewise(Expr* res);

void cond_init(void);

#endif

#ifndef MATHILDA_TRIGSIMP_H
#define MATHILDA_TRIGSIMP_H

#include "expr.h"

Expr* builtin_trigtoexp(Expr* res);
Expr* builtin_exptotrig(Expr* res);
Expr* builtin_trigexpand(Expr* res);
Expr* builtin_trigfactor(Expr* res);
Expr* builtin_trigreduce(Expr* res);

void trigsimp_init(void);

#endif /* MATHILDA_TRIGSIMP_H */

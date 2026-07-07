#ifndef JOIN_H
#define JOIN_H

#include "expr.h"

Expr* builtin_join(Expr* res);
Expr* builtin_catenate(Expr* res);

#endif /* JOIN_H */

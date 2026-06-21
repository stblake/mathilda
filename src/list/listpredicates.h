#ifndef LISTPREDICATES_H
#define LISTPREDICATES_H

#include "expr.h"

Expr* builtin_listq(Expr* res);
Expr* builtin_vectorq(Expr* res);
Expr* builtin_matrixq(Expr* res);

#endif /* LISTPREDICATES_H */

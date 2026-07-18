#ifndef SETOPS_H
#define SETOPS_H

#include "expr.h"

Expr* builtin_tally(Expr* res);
Expr* builtin_union(Expr* res);
Expr* builtin_intersection(Expr* res);
Expr* builtin_deleteduplicates(Expr* res);
Expr* builtin_deleteduplicatesby(Expr* res);
Expr* builtin_commonest(Expr* res);

#endif /* SETOPS_H */

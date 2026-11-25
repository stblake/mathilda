#ifndef DEFAULT_HELPER_H
#define DEFAULT_HELPER_H

#include "expr.h"

Expr* get_default_value(Expr* pat_head, int pos, int total_pats);
/* Default[f] / Default[f, i] for the built-in heads (Plus, Times, Power). */
Expr* builtin_default(Expr* res);

#endif // DEFAULT_HELPER_H

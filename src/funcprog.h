
#ifndef FUNCPROG_H
#define FUNCPROG_H

#include "expr.h"

Expr* builtin_apply(Expr* res);
Expr* builtin_map(Expr* res);
Expr* builtin_mapindexed(Expr* res);
Expr* builtin_map_all(Expr* res);
Expr* builtin_map_at(Expr* res);
Expr* builtin_nest(Expr* res);
Expr* builtin_nestlist(Expr* res);
Expr* builtin_nestwhile(Expr* res);
Expr* builtin_nestwhilelist(Expr* res);
Expr* builtin_fixedpointlist(Expr* res);
Expr* builtin_fixedpoint(Expr* res);
Expr* builtin_fold(Expr* res);
Expr* builtin_foldlist(Expr* res);
Expr* builtin_select(Expr* res);
Expr* builtin_select_first(Expr* res);
Expr* builtin_takewhile(Expr* res);
Expr* builtin_lengthwhile(Expr* res);
Expr* builtin_scan(Expr* res);
Expr* builtin_catch(Expr* res);
Expr* builtin_throw(Expr* res);
Expr* builtin_all_true(Expr* res);
Expr* builtin_any_true(Expr* res);
Expr* builtin_none_true(Expr* res);
Expr* builtin_through(Expr* res);
Expr* builtin_thread(Expr* res);
Expr* builtin_freeq(Expr* res);
Expr* builtin_distribute(Expr* res);

#endif

Expr* builtin_inner(Expr* res);
Expr* builtin_outer(Expr* res);
Expr* builtin_tuples(Expr* res);
Expr* builtin_permutations(Expr* res);

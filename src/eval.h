#ifndef EVAL_H
#define EVAL_H

#include "expr.h"
#include "attr.h"
#include <stdint.h>

// Single step evaluation
Expr* evaluate_step(Expr* e);

// Full infinite evaluation loop (stops at fixed point)
Expr* evaluate(Expr* e);

// Recursion-depth limit for nested evaluate() calls. Guards C-stack overflow
// when expressions trigger deeply recursive sub-evaluation (Listable, Flat,
// user DownValues that re-evaluate, etc.). Default 1024 (matches Mathematica).
int  eval_get_recursion_limit(void);
void eval_set_recursion_limit(int n);
int  eval_get_recursion_depth(void);

// Seeds the user-visible $RecursionLimit symbol with its default value
// and installs the docstring. Call after symtab_init().
void eval_init(void);

// Helper to evaluate and free the input expression
static inline Expr* eval_and_free(Expr* e) {
    if (!e) return NULL;
    Expr* res = evaluate(e);
    expr_free(e);
    return res;
}

// Attributes built-in
Expr* builtin_attributes(Expr* res);

void eval_flatten_args(Expr* e, const char* head_name);
int eval_compare_expr_ptrs(const void* a, const void* b);
static inline void eval_sort_args(Expr* e) { if (e->type == EXPR_FUNCTION && e->data.function.arg_count > 0) qsort(e->data.function.args, e->data.function.arg_count, sizeof(Expr*), eval_compare_expr_ptrs); }
#endif // EVAL_H

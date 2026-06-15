/* sum_util.c -- shared helpers for the Sum sub-algorithms. */

#include "sum_internal.h"
#include "sym_names.h"
#include "eval.h"
#include "expr.h"
#include <string.h>
#include <stdlib.h>

Expr* sum_eval(const char* head, Expr** args, size_t n) {
    Expr* call = expr_new_function(expr_new_symbol(head), args, n);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

Expr* sum_int(int64_t v) { return expr_new_integer(v); }

Expr* sum_subst(Expr* e, Expr* var, Expr* val) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                    (Expr*[]){ expr_copy(var), expr_copy(val) }, 2);
    Expr* args[2] = { expr_copy(e), rule };
    return sum_eval("ReplaceAll", args, 2);
}

Expr* sum_factor(Expr* e) {
    Expr* arg[1] = { expr_copy(e) };
    Expr* r = sum_eval("Factor", arg, 1);
    /* If Factor came back unevaluated (head Factor), keep the input form. */
    if (r && r->type == EXPR_FUNCTION
          && r->data.function.head->type == EXPR_SYMBOL
          && strcmp(r->data.function.head->data.symbol, "Factor") == 0) {
        expr_free(r);
        return expr_copy(e);
    }
    return r;
}

bool sum_free_of(Expr* e, Expr* var) {
    Expr* args[2] = { expr_copy(e), expr_copy(var) };
    Expr* r = sum_eval("FreeQ", args, 2);
    bool yes = (r && r->type == EXPR_SYMBOL && strcmp(r->data.symbol, "True") == 0);
    if (r) expr_free(r);
    return yes;
}

Expr* sum_sub(Expr* a, Expr* b) {
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(b) }, 2);
    Expr* args[2] = { expr_copy(a), neg };
    return sum_eval("Plus", args, 2);
}

bool sum_stage_args(Expr* res, Expr** f, Expr** var,
                    Expr** imin, Expr** imax, bool* definite) {
    if (res->type != EXPR_FUNCTION) return false;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;
    if (argc != 2 && argc != 4) return false;
    if (a[1]->type != EXPR_SYMBOL) return false;   /* summation var */
    *f = a[0];
    *var = a[1];
    if (argc == 4) {
        *imin = a[2];
        *imax = a[3];
        *definite = true;
    } else {
        *imin = NULL;
        *imax = NULL;
        *definite = false;
    }
    return true;
}

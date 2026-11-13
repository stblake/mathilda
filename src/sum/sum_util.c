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
          && r->data.function.head->data.symbol.name == SYM_Factor) {
        expr_free(r);
        return expr_copy(e);
    }
    return r;
}

bool sum_free_of(Expr* e, Expr* var) {
    Expr* args[2] = { expr_copy(e), expr_copy(var) };
    Expr* r = sum_eval("FreeQ", args, 2);
    bool yes = (r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True);
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

/* Head test on an interned name. */
static bool su_head(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}

/* Record the integer root of one irreducible factor b of a denominator, when
 * b is linear in var with exact numeric coefficients: b = a var + c gives
 * var = -c/a. */
static void su_linear_integer_root(Expr* b, Expr* var, int64_t* out,
                                   size_t cap, size_t* n) {
    if (*n >= cap || sum_free_of(b, var)) return;
    Expr* deg = sum_eval("Exponent", (Expr*[]){ expr_copy(b), expr_copy(var) }, 2);
    bool linear = deg && deg->type == EXPR_INTEGER && deg->data.integer == 1;
    if (deg) expr_free(deg);
    if (!linear) return;
    Expr* a = sum_eval("Coefficient", (Expr*[]){ expr_copy(b), expr_copy(var), sum_int(1) }, 3);
    Expr* c = sum_eval("Coefficient", (Expr*[]){ expr_copy(b), expr_copy(var), sum_int(0) }, 3);
    Expr* r = sum_eval("Times", (Expr*[]){ sum_int(-1), c,
                  expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ a, sum_int(-1) }, 2) }, 3);
    if (r && r->type == EXPR_INTEGER) out[(*n)++] = r->data.integer;
    if (r) expr_free(r);
}

size_t sum_integer_poles(Expr* f, Expr* var, int64_t* out, size_t cap) {
    Expr* tg = sum_eval("Together", (Expr*[]){ expr_copy(f) }, 1);
    Expr* den = sum_eval("Denominator", (Expr*[]){ tg }, 1);
    size_t n = 0;
    if (!den || sum_free_of(den, var)) { if (den) expr_free(den); return 0; }
    Expr* fac = sum_factor(den);
    expr_free(den);
    size_t nf = su_head(fac, SYM_Times) ? fac->data.function.arg_count : 1;
    for (size_t i = 0; i < nf; i++) {
        Expr* t = su_head(fac, SYM_Times) ? fac->data.function.args[i] : fac;
        if (su_head(t, SYM_Power) && t->data.function.arg_count == 2
            && t->data.function.args[1]->type == EXPR_INTEGER
            && t->data.function.args[1]->data.integer > 0)
            t = t->data.function.args[0];
        su_linear_integer_root(t, var, out, cap, &n);
    }
    expr_free(fac);
    return n;
}

/*
 * product_special.c -- Product`Special (Stage 6): named special-function
 * product identities.
 *
 *   prod_{k=1}^{n}   k^k     = Hyperfactorial[n],
 *   prod_{k=1}^{n-1} Gamma[k] = BarnesG[n]   (i.e. prod_{k=1}^{m} Gamma[k]
 *                                             = BarnesG[m+1]).
 *
 * Recognised structurally in C (the held Plus/Power forms defeat .m DownValue
 * matching for these shapes), only at lower bound 1.
 */

#include "product_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "sym_names.h"
#include <string.h>

static bool is_one(const Expr* e) {
    return e->type == EXPR_INTEGER && e->data.integer == 1;
}

/* Power[var, var]  (i.e. k^k). */
static bool is_k_pow_k(const Expr* f, const Expr* var) {
    return f->type == EXPR_FUNCTION
        && f->data.function.head->type == EXPR_SYMBOL
        && f->data.function.head->data.symbol == SYM_Power
        && f->data.function.arg_count == 2
        && f->data.function.args[0]->type == EXPR_SYMBOL
        && f->data.function.args[0]->data.symbol == var->data.symbol
        && f->data.function.args[1]->type == EXPR_SYMBOL
        && f->data.function.args[1]->data.symbol == var->data.symbol;
}

/* Gamma[var]. */
static bool is_gamma_k(const Expr* f, const Expr* var) {
    return f->type == EXPR_FUNCTION
        && f->data.function.head->type == EXPR_SYMBOL
        && strcmp(f->data.function.head->data.symbol, "Gamma") == 0
        && f->data.function.arg_count == 1
        && f->data.function.args[0]->type == EXPR_SYMBOL
        && f->data.function.args[0]->data.symbol == var->data.symbol;
}

Expr* builtin_product_special(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!product_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite || !is_one(imin)) return NULL;
    if (imax->type == EXPR_SYMBOL && imax->data.symbol == SYM_Infinity) return NULL;

    /* prod_{k=1}^{n} k^k = Hyperfactorial[n]. */
    if (is_k_pow_k(f, var)) {
        return expr_new_function(expr_new_symbol("Hyperfactorial"),
                   (Expr*[]){ expr_copy(imax) }, 1);
    }

    /* prod_{k=1}^{n} Gamma[k] = BarnesG[n+1]. */
    if (is_gamma_k(f, var)) {
        Expr* np1 = expr_new_function(expr_new_symbol(SYM_Plus),
                        (Expr*[]){ expr_copy(imax), expr_new_integer(1) }, 2);
        Expr* arg = evaluate(np1);
        expr_free(np1);
        return expr_new_function(expr_new_symbol("BarnesG"), &arg, 1);  /* adopts arg */
    }

    return NULL;
}

void product_special_init(void) {
    symtab_add_builtin("Product`Special", builtin_product_special);
    symtab_get_def("Product`Special")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Product`Special",
        "Product`Special[f, i, imin, imax] recognises named special-function "
        "products: prod i^i = Hyperfactorial, prod Gamma[i] = BarnesG. Returns "
        "unevaluated otherwise.");
}

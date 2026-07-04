/*
 * product_besselzero.c -- Product`BesselZero: Hadamard product over Bessel zeros.
 *
 * The entire-function Hadamard product for the Bessel function of the first
 * kind (DLMF 10.21.15) is
 *
 *   J_nu(x) = (x/2)^nu / Gamma[nu+1] * prod_{k=1}^Inf (1 - x^2 / j_{nu,k}^2),
 *
 * where j_{nu,k} = BesselJZero[nu, k].  Hence
 *
 *   Product[1 - x^2/BesselJZero[nu,k]^2, {k,1,Inf}]
 *       = Gamma[nu+1] (2/x)^nu BesselJ[nu, x],
 *
 * which for nu = 0 is simply BesselJ[0, x].
 *
 * Recognises a body  1 - X2 * BesselJZero[nu, k]^(-2)  with X2 (= x^2) and nu
 * free of k, the lower bound imin == 1 (the canonical Hadamard lower bound),
 * and emits the closed form above (letting evaluate() collapse the nu==0
 * prefactor Gamma[1] (2/x)^0 = 1).
 *
 * Memory contract: takes ownership of res but must not free it; returns an
 * owned closed form or NULL to fall through.
 */

#include "product_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>

static bool bz_is_inf(const Expr* e) {
    return e->type == EXPR_SYMBOL && e->data.symbol == SYM_Infinity;
}
static bool bz_is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, name) == 0;
}

/* Is e == Power[BesselJZero[nu, var], -2] with nu free of var?  Sets *nu_out
 * (alias) on success. */
static bool bz_match_jzero_pow(Expr* e, Expr* var, Expr** nu_out) {
    if (!bz_is_head(e, "Power") || e->data.function.arg_count != 2) return false;
    Expr* base = e->data.function.args[0];
    Expr* exp  = e->data.function.args[1];
    if (!(exp->type == EXPR_INTEGER && exp->data.integer == -2)) return false;
    if (!bz_is_head(base, "BesselJZero") || base->data.function.arg_count != 2) return false;
    Expr* nu = base->data.function.args[0];
    Expr* kk = base->data.function.args[1];
    if (!expr_eq(kk, var)) return false;
    if (!prod_free_of(nu, var)) return false;
    *nu_out = nu;
    return true;
}

Expr* builtin_product_besselzero(Expr* res);

Expr* builtin_product_besselzero(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!product_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite || !bz_is_inf(imax)) return NULL;
    if (!(imin->type == EXPR_INTEGER && imin->data.integer == 1)) return NULL;

    Expr* fn = evaluate(expr_copy(f));

    /* f must be Plus[1, t] with t = -X2 * BesselJZero[nu,k]^-2. */
    if (!bz_is_head(fn, "Plus") || fn->data.function.arg_count != 2) {
        expr_free(fn); return NULL;
    }
    /* Identify the integer-1 term and the other term t. */
    Expr *one_term = NULL, *t = NULL;
    for (int i = 0; i < 2; i++) {
        Expr* a = fn->data.function.args[i];
        if (a->type == EXPR_INTEGER && a->data.integer == 1) one_term = a;
        else t = a;
    }
    if (!one_term || !t) { expr_free(fn); return NULL; }

    /* sub = -t = X2 * BesselJZero[nu,k]^-2. */
    Expr* subraw = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), expr_copy(t) }, 2);
    Expr* sub = evaluate(subraw);
    expr_free(subraw);

    /* Scan sub's factors: exactly one BesselJZero[nu,k]^-2, the rest = X2 (free of k). */
    bool is_times = bz_is_head(sub, "Times");
    size_t nf = is_times ? sub->data.function.arg_count : 1;
    Expr* nu = NULL;
    Expr* x2 = expr_new_integer(1);
    bool ok = true, found = false;
    for (size_t i = 0; i < nf && ok; i++) {
        Expr* fac = is_times ? sub->data.function.args[i] : sub;
        Expr* nutmp = NULL;
        if (bz_match_jzero_pow(fac, var, &nutmp)) {
            if (found) { ok = false; break; }
            found = true; nu = nutmp;
        } else if (prod_free_of(fac, var)) {
            Expr* p = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ x2, expr_copy(fac) }, 2);   /* adopts x2 */
            x2 = evaluate(p); expr_free(p);
        } else { ok = false; }
    }

    if (!ok || !found) {
        expr_free(fn); expr_free(sub); expr_free(x2); return NULL;
    }

    /* xroot = PowerExpand[Sqrt[X2]] (= x for X2 = x^2). */
    Expr* sq = expr_new_function(expr_new_symbol("Sqrt"),
                   (Expr*[]){ expr_copy(x2) }, 1);
    Expr* xroot = prod_eval("PowerExpand", (Expr*[]){ sq }, 1);   /* adopts sq */
    Expr* nu_copy = expr_copy(nu);
    expr_free(x2);
    expr_free(fn); expr_free(sub);   /* nu aliased into sub; copied above */

    /* Gamma[nu+1] * (2/xroot)^nu * BesselJ[nu, xroot]. */
    Expr* nup1 = expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ expr_copy(nu_copy), expr_new_integer(1) }, 2);
    Expr* gam = expr_new_function(expr_new_symbol("Gamma"),
                    (Expr*[]){ nup1 }, 1);
    Expr* two = expr_new_integer(2);
    Expr* two_over_x = prod_div(two, xroot);   /* prod_div copies both args */
    expr_free(two);
    Expr* pref = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ two_over_x, expr_copy(nu_copy) }, 2);
    Expr* bj = expr_new_function(expr_new_symbol(SYM_BesselJ),
                   (Expr*[]){ nu_copy, expr_copy(xroot) }, 2);   /* adopts nu_copy */
    Expr* full = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ gam, pref, bj }, 3);
    expr_free(xroot);
    Expr* out = evaluate(full);
    expr_free(full);
    return out;
}

void product_besselzero_init(void) {
    symtab_add_builtin("Product`BesselZero", builtin_product_besselzero);
    symtab_get_def("Product`BesselZero")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Product`BesselZero",
        "Product`BesselZero[f, i, 1, Infinity] evaluates the Bessel Hadamard "
        "product Product[1 - x^2/BesselJZero[n,i]^2] = "
        "Gamma[n+1] (2/x)^n BesselJ[n, x].");
}

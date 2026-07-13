/*
 * product_logsum.c -- Product`LogSum: the general Exp/log-sum bridge.
 *
 * The universal identity prod_k f(k) = Exp[ Sum_k Log[f(k)] ] closes an
 * infinite product whenever Log[f(k)] decomposes into pieces the shipped Sum
 * can sum:
 *
 *   Product[Exp[(-1)^k/k], {k,1,Inf}] = Exp[Sum[(-1)^k/k]] = Exp[-Log 2] = 1/2,
 *   Product[k^(1/k^2),     {k,1,Inf}] = Exp[Sum[Log[k]/k^2]] = Exp[-Zeta'[2]],
 *   Product[E^(1/k)(1-1/k),{k,2,Inf}] = Exp[Sum[1/k+Log[(k-1)/k]]] = E^(g-1).
 *
 * PowerExpand turns Log[Exp[g]] -> g and Log of a product/power into a sum of
 * logs (plain Log[Exp[.]] does NOT auto-reduce, so PowerExpand is required).
 *
 * This is the most general and most expensive route, so it runs LAST in the
 * cascade: it only ever sees products no structural stage could close.  It
 * engages only when the body carries a symbolic power (Exp[.] is Power[E,.]),
 * which excludes plain rationals (handled earlier) -- taking logs of a rational
 * and calling Sum would churn for nothing.
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

static bool is_inf_sym(const Expr* e) {
    return e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity;
}

static bool head_sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

Expr* builtin_product_logsum(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!product_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite || !is_inf_sym(imax)) return NULL;

    /* Engage only on a genuine symbolic-power body (Exp[.] == Power[E,.]); pure
     * rationals were handled by earlier stages. */
    if (!prod_has_symbolic_power(f, var)) return NULL;

    /* L = PowerExpand[Log[f]]. */
    Expr* logf = expr_new_function(expr_new_symbol("Log"),
                     (Expr*[]){ expr_copy(f) }, 1);
    Expr* L = prod_eval("PowerExpand", (Expr*[]){ logf }, 1);   /* adopts logf */

    /* If Log did not decompose (L is still literally Log[f]), give up -- the
     * Sum would just churn. */
    Expr* logf_ref = expr_new_function(expr_new_symbol("Log"),
                         (Expr*[]){ expr_copy(f) }, 1);
    bool undecomposed = expr_eq(L, logf_ref);
    expr_free(logf_ref);
    if (undecomposed) { expr_free(L); return NULL; }

    /* S = Sum[L, {var, imin, Infinity}]. */
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List),
                     (Expr*[]){ expr_copy(var), expr_copy(imin),
                                expr_copy(imax) }, 3);
    Expr* S = prod_eval("Sum", (Expr*[]){ L, spec }, 2);        /* adopts L, spec */

    /* Sum must have closed (head no longer Sum) and be free of the index. */
    if (head_sym_is(S, "Sum") || !prod_free_of(S, var)) { expr_free(S); return NULL; }

    /* Result = Exp[S]. */
    Expr* out = expr_new_function(expr_new_symbol("Exp"), (Expr*[]){ S }, 1); /* adopts S */
    Expr* r = evaluate(out);
    expr_free(out);
    return r;
}

void product_logsum_init(void) {
    symtab_add_builtin("Product`LogSum", builtin_product_logsum);
    symtab_get_def("Product`LogSum")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Product`LogSum",
        "Product`LogSum[f, i, imin, Infinity] evaluates prod f(i) via "
        "Exp[Sum[PowerExpand[Log[f(i)]], {i,imin,Infinity}]] when the log-sum "
        "closes. Engages only on symbolic-power bodies; returns unevaluated "
        "otherwise.");
}

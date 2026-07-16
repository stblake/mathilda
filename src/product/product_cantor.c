/*
 * product_cantor.c -- Product`Cantor: double-exponential (Cantor) telescoping.
 *
 * The identity (1 - x) prod_{k=0}^Inf (1 + x^(2^k)) = 1 gives, for |x| < 1,
 *
 *   prod_{k=imin}^Inf (1 + x^(2^L(k))) = 1 / (1 - x^(2^L(imin)))
 *
 * whenever the double-exponent doubles each step (2^L(k+1) == 2 * 2^L(k), i.e.
 * L is affine in k).  Concretely:
 *
 *   Product[1 + (1/3)^(2^k), {k, 0, Infinity}]  ->  1/(1 - 1/3) = 3/2.
 *
 * Gate (all required): infinite upper bound; body == 1 + base^E with base free
 * of k; E involves k and doubles each step (Simplify[E(k+1)/E(k)] == 2);
 * Abs[base] < 1 is provably True.  Anything else -> NULL.
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

static bool is_power(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2;
}

static bool is_int_val(const Expr* e, int64_t v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}

/* Evaluate the boolean predicate Less[Abs[base], 1] and test for literal True. */
static bool abs_lt_one(Expr* base) {
    Expr* ab = expr_new_function(expr_new_symbol("Abs"),
                   (Expr*[]){ expr_copy(base) }, 1);
    Expr* lt = expr_new_function(expr_new_symbol(SYM_Less),
                   (Expr*[]){ ab, expr_new_integer(1) }, 2);   /* adopts ab */
    Expr* r = evaluate(lt);
    expr_free(lt);
    bool yes = (r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True);
    if (r) expr_free(r);
    return yes;
}

Expr* builtin_product_cantor(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!product_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite || !is_inf_sym(imax)) return NULL;

    /* Body must be Plus with exactly two args, one of which is the literal 1. */
    if (!(f->type == EXPR_FUNCTION
          && f->data.function.head->type == EXPR_SYMBOL
          && f->data.function.head->data.symbol.name == SYM_Plus
          && f->data.function.arg_count == 2)) return NULL;

    Expr* a0 = f->data.function.args[0];
    Expr* a1 = f->data.function.args[1];
    Expr* term;
    if (is_int_val(a0, 1))      term = a1;
    else if (is_int_val(a1, 1)) term = a0;
    else return NULL;

    /* term must be base^E, base free of k, E involving k. */
    if (!is_power(term)) return NULL;
    Expr* base = term->data.function.args[0];
    Expr* E    = term->data.function.args[1];
    if (!prod_free_of(base, var)) return NULL;
    if (prod_free_of(E, var))     return NULL;

    /* Convergence: |base| < 1 must be provably true. */
    if (!abs_lt_one(base)) return NULL;

    /* Double-exponent doubling test: Simplify[E(k+1)/E(k)] == 2. */
    Expr* kp1 = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(var), expr_new_integer(1) }, 2);
    Expr* E_next = prod_subst(E, var, kp1);
    expr_free(kp1);
    Expr* ratio_raw = prod_div(E_next, E);
    expr_free(E_next);
    Expr* ratio = prod_eval("Simplify", (Expr*[]){ ratio_raw }, 1);  /* adopts */
    bool doubling = is_int_val(ratio, 2);
    expr_free(ratio);
    if (!doubling) return NULL;

    /* Closed form: M = E(imin);  result = 1 / (1 - base^M). */
    Expr* M = prod_subst(E, var, imin);
    Expr* baseM = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_copy(base), M }, 2);   /* adopts M */
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), baseM }, 2);  /* adopts baseM */
    Expr* denom = expr_new_function(expr_new_symbol(SYM_Plus),
                      (Expr*[]){ expr_new_integer(1), neg }, 2);    /* adopts neg */
    Expr* one = expr_new_integer(1);
    Expr* out = prod_div(one, denom);
    expr_free(one);
    expr_free(denom);
    return out;
}

void product_cantor_init(void) {
    symtab_add_builtin("Product`Cantor", builtin_product_cantor);
    symtab_get_def("Product`Cantor")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Product`Cantor",
        "Product`Cantor[f, i, imin, Infinity] evaluates double-exponential "
        "(Cantor) products prod (1 + x^(2^i)) with |x|<1 via 1/(1 - x^(2^i(imin))). "
        "Returns unevaluated otherwise.");
}

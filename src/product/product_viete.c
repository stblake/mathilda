/*
 * product_viete.c -- Product`Viete: Viete-type cosine products via the
 * double-angle telescoping identity.
 *
 * For a body Cos[a(k)] whose angle halves each step (a(k+1) == a(k)/2), the
 * infinite product telescopes:
 *
 *   prod_{k=imin}^Inf Cos[a(k)]
 *     = prod_{j=0}^Inf Cos[y / 2^j]              (y = a(imin))
 *     = Sin[2 y] / (2 y),
 *
 * because prod_{j=1}^N cos(x/2^j) = sin(x)/(2^N sin(x/2^N)) -> sin(x)/x, and
 * shifting the index so the product starts at j=0 replaces x by 2y.
 *
 *   Product[Cos[Pi/2^(k+1)], {k, 1, Infinity}]  ->  Sin[Pi/2]/(Pi/2) = 2/Pi.
 *
 * The angle-halving test is representation-independent: Simplify[a(k+1)/a(k)]
 * must be exactly 1/2.  This also guarantees convergence (the factors -> 1).
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

/* True iff e is exactly the rational 1/2 (Rational[1,2]). */
static bool is_one_half(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_INTEGER
        && e->data.function.args[0]->data.integer == 1
        && e->data.function.args[1]->type == EXPR_INTEGER
        && e->data.function.args[1]->data.integer == 2;
}

Expr* builtin_product_viete(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!product_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite || !is_inf_sym(imax)) return NULL;

    /* Body must be a single Cos[angle] that genuinely depends on k. */
    if (!head_sym_is(f, "Cos") || f->data.function.arg_count != 1) return NULL;
    Expr* angle = f->data.function.args[0];
    if (prod_free_of(angle, var)) return NULL;

    /* Angle-halving test: Simplify[angle(k+1)/angle(k)] == 1/2. */
    Expr* kp1 = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(var), expr_new_integer(1) }, 2);
    Expr* a_next = prod_subst(angle, var, kp1);
    expr_free(kp1);
    Expr* ratio_raw = prod_div(a_next, angle);      /* a(k+1) / a(k) */
    expr_free(a_next);
    Expr* ratio = prod_eval("Simplify", (Expr*[]){ ratio_raw }, 1);  /* adopts */
    bool halving = is_one_half(ratio);
    expr_free(ratio);
    if (!halving) return NULL;

    /* Closed form: y = angle(imin), x = 2 y, result = Sin[x]/x. */
    Expr* y = prod_subst(angle, var, imin);
    Expr* x = expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ expr_new_integer(2), y }, 2);   /* adopts y */
    Expr* xe = evaluate(x);
    expr_free(x);
    Expr* sinx = expr_new_function(expr_new_symbol("Sin"),
                     (Expr*[]){ expr_copy(xe) }, 1);
    Expr* out = prod_div(sinx, xe);
    expr_free(sinx);
    expr_free(xe);
    return out;
}

void product_viete_init(void) {
    symtab_add_builtin("Product`Viete", builtin_product_viete);
    symtab_get_def("Product`Viete")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Product`Viete",
        "Product`Viete[f, i, imin, Infinity] evaluates Viete-type cosine "
        "products prod Cos[a(i)] whose angle halves each step "
        "(a(i+1)==a(i)/2) via Sin[2 a(imin)]/(2 a(imin)). Returns unevaluated "
        "otherwise.");
}

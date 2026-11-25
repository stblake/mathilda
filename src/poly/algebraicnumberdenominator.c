/* algebraicnumberdenominator.c — AlgebraicNumberDenominator[a].
 *
 * See algebraicnumberdenominator.h. The value — the smallest positive integer n
 * with n a an algebraic integer — is computed exactly by
 * flint_qqbar_algebraic_number_denominator (a per-prime valuation over the
 * minimal polynomial). This file only checks the argument shape, maps the
 * tri-state engine result to Integer / message / unevaluated, and (being
 * Listable) lets the evaluator thread over a list of algebraic numbers.
 */

#include "algebraicnumberdenominator.h"

#include "flint_qqbar.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "eval.h"
#include "sym_names.h"

#include <stdio.h>
#include <stdlib.h>

/* Collect the denominators of the Rational atoms of e into *acc (an LCM
 * call's argument list). */
static void and_collect_dens(const Expr* e, Expr*** acc, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational &&
        e->data.function.arg_count == 2) {
        if (*n == *cap) { *cap = *cap ? 2 * *cap : 8; *acc = realloc(*acc, *cap * sizeof(Expr*)); }
        (*acc)[(*n)++] = expr_copy(e->data.function.args[1]);
        return;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        and_collect_dens(e->data.function.args[i], acc, n, cap);
}

/* Root[f, k] with rational coefficients as Root[L f, k] with integer ones
 * (L the LCM of the denominators): the same number, in the form the engine
 * reads (audit AlgebraicNumberDenominator#1).  NULL when not applicable. */
static Expr* and_integer_root(const Expr* a) {
    if (a->type != EXPR_FUNCTION || a->data.function.head->type != EXPR_SYMBOL ||
        a->data.function.head->data.symbol.name != SYM_Root || a->data.function.arg_count != 2)
        return NULL;
    const Expr* fn = a->data.function.args[0];
    if (fn->type != EXPR_FUNCTION || fn->data.function.head->type != EXPR_SYMBOL ||
        fn->data.function.head->data.symbol.name != SYM_Function || fn->data.function.arg_count != 1)
        return NULL;
    Expr** dens = NULL; size_t nd = 0, cap = 0;
    and_collect_dens(fn->data.function.args[0], &dens, &nd, &cap);
    if (nd == 0) { free(dens); return NULL; }
    Expr* L = eval_and_free(expr_new_function(expr_new_symbol(SYM_LCM), dens, nd));
    free(dens);
    Expr* body = eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand),
                     (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Times),
                         (Expr*[]){ L, expr_copy(fn->data.function.args[0]) }, 2) }, 1));
    Expr* nf = expr_new_function(expr_new_symbol(SYM_Function), &body, 1);
    return expr_new_function(expr_new_symbol(SYM_Root),
                             (Expr*[]){ nf, expr_copy(a->data.function.args[1]) }, 2);
}

Expr* builtin_algebraicnumberdenominator(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;

    const Expr* a = res->data.function.args[0];
    Expr* out = NULL;
    int r = flint_qqbar_algebraic_number_denominator(a, &out);
    if (r == 0) {
        Expr* ir = and_integer_root(a);
        if (ir) {
            Expr* out2 = NULL;
            int r2 = flint_qqbar_algebraic_number_denominator(ir, &out2);
            expr_free(ir);
            if (r2 > 0) return out2;
        }
    }
    if (r > 0) return out;                   /* success: the denominator */
    if (r < 0) return NULL;                  /* undecided (FLINT compiled out) */

    /* r == 0: not a constant algebraic number — report and stay unevaluated. */
    char* s = expr_to_string((Expr*)a);
    mth_message("AlgebraicNumberDenominator", "nalg", "%s is not an explicit algebraic number.", s ? s : "?");
    free(s);
    return NULL;
}

void algebraicnumberdenominator_init(void) {
    symtab_add_builtin("AlgebraicNumberDenominator",
                       builtin_algebraicnumberdenominator);
    /* Listable: AlgebraicNumberDenominator[{a1, ..}] threads element-wise, giving
     * one denominator per algebraic number. */
    symtab_get_def("AlgebraicNumberDenominator")->attributes |=
        ATTR_LISTABLE | ATTR_PROTECTED;
}

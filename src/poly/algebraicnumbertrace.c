/* algebraicnumbertrace.c — AlgebraicNumberTrace[a].
 *
 * See algebraicnumbertrace.h.  The value — the field trace of an algebraic number,
 * optionally relative to a field Q(theta) given by Extension -> theta — is
 * computed exactly by flint_qqbar_algebraic_number_trace (both cases reduce to a
 * minimal-polynomial read; see flint_qqbar.c).  This file only separates the
 * positional argument from the Extension option, checks arity, maps the
 * tri/quad-state engine result to Integer/Rational / message / unevaluated, and
 * (being Listable) lets the evaluator thread over a list of algebraic numbers.
 */

#include "algebraicnumbertrace.h"

#include "flint_qqbar.h"
#include "options.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "sym_names.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */

#include <stdio.h>
#include <stdlib.h>

Expr* builtin_algebraicnumbertrace(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    /* Extension -> theta is a trailing option; NULL means None/Automatic/absent,
     * i.e. the absolute trace over Q.  argc is left holding the positional count. */
    size_t argc = res->data.function.arg_count;
    const Expr* theta = extract_extension_option(res, &argc);
    if (argc != 1) return NULL;
    /* AlgebraicNumberTrace#2: of repeated Extension options the first is the
     * one in force, as for every option (the shared reader keeps the last). */
    for (size_t i = 1; i < res->data.function.arg_count; i++) {
        const Expr* o = res->data.function.args[i];
        if (o->type == EXPR_FUNCTION && o->data.function.arg_count == 2 &&
            o->data.function.head->type == EXPR_SYMBOL &&
            (o->data.function.head->data.symbol.name == SYM_Rule ||
             o->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
            o->data.function.args[0]->type == EXPR_SYMBOL &&
            o->data.function.args[0]->data.symbol.name == SYM_Extension) {
            const Expr* v = o->data.function.args[1];
            bool trivial = v->type == EXPR_SYMBOL &&
                           (v->data.symbol.name == SYM_Automatic || v->data.symbol.name == SYM_None);
            theta = trivial ? NULL : v;
            break;
        }
    }

    const Expr* a = res->data.function.args[0];
    Expr* out = NULL;
    int r = flint_qqbar_algebraic_number_trace(a, theta, &out);
    if (r == 1) return out;                   /* success: the trace */
    if (r < 0) return NULL;                    /* undecided (FLINT compiled out) */

    if (r == 2) {                              /* a is not an element of Q(theta) */
        char* s = expr_to_string((Expr*)a);
        mth_message("AlgebraicNumberTrace", "ext", "%s is not an element of the field extension specified by the Extension option.", s ? s : "?");
        free(s);
        return NULL;
    }

    /* r == 0: not a constant algebraic number — report and stay unevaluated.
     * When a is one, the culprit is the Extension value (AlgebraicNumberNorm#2,
     * AlgebraicNumberTrace#1: Extension -> Pi blamed Sqrt[2]). */
    const Expr* culprit = (theta && flint_qqbar_is_constant_algebraic(a)) ? theta : a;
    char* s = expr_to_string((Expr*)culprit);
    mth_message("AlgebraicNumberTrace", "nalg", "%s is not an explicit algebraic number.", s ? s : "?");
    free(s);
    return NULL;
}

void algebraicnumbertrace_init(void) {
    symtab_add_builtin("AlgebraicNumberTrace", builtin_algebraicnumbertrace);
    /* Listable: AlgebraicNumberTrace[{a1, ..}] threads element-wise, giving one
     * trace per algebraic number; a trailing Extension option is repeated. */
    symtab_get_def("AlgebraicNumberTrace")->attributes |=
        ATTR_LISTABLE | ATTR_PROTECTED;
}

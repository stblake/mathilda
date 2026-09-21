/* algebraicnumbernorm.c — AlgebraicNumberNorm[a].
 *
 * See algebraicnumbernorm.h.  The value — the field norm of an algebraic number,
 * optionally relative to a field Q(theta) given by Extension -> theta — is
 * computed exactly by flint_qqbar_algebraic_number_norm (both cases reduce to a
 * minimal-polynomial read; see flint_qqbar.c).  This file only separates the
 * positional argument from the Extension option, checks arity, maps the
 * tri/quad-state engine result to Integer/Rational / message / unevaluated, and
 * (being Listable) lets the evaluator thread over a list of algebraic numbers.
 */

#include "algebraicnumbernorm.h"

#include "flint_qqbar.h"
#include "options.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>

Expr* builtin_algebraicnumbernorm(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    /* Extension -> theta is a trailing option; NULL means None/Automatic/absent,
     * i.e. the absolute norm over Q.  argc is left holding the positional count. */
    size_t argc = res->data.function.arg_count;
    const Expr* theta = extract_extension_option(res, &argc);
    if (argc != 1) return NULL;

    const Expr* a = res->data.function.args[0];
    Expr* out = NULL;
    int r = flint_qqbar_algebraic_number_norm(a, theta, &out);
    if (r == 1) return out;                  /* success: the norm */
    if (r < 0) return NULL;                   /* undecided (FLINT compiled out) */

    if (r == 2) {                             /* a is not an element of Q(theta) */
        char* s = expr_to_string((Expr*)a);
        fprintf(stderr,
                "AlgebraicNumberNorm::ext: %s is not an element of the field "
                "extension specified by the Extension option.\n", s ? s : "?");
        free(s);
        return NULL;
    }

    /* r == 0: not a constant algebraic number — report and stay unevaluated. */
    char* s = expr_to_string((Expr*)a);
    fprintf(stderr,
            "AlgebraicNumberNorm::nalg: %s is not an explicit algebraic number.\n",
            s ? s : "?");
    free(s);
    return NULL;
}

void algebraicnumbernorm_init(void) {
    symtab_add_builtin("AlgebraicNumberNorm", builtin_algebraicnumbernorm);
    /* Listable: AlgebraicNumberNorm[{a1, ..}] threads element-wise, giving one
     * norm per algebraic number; a trailing Extension option is repeated. */
    symtab_get_def("AlgebraicNumberNorm")->attributes |=
        ATTR_LISTABLE | ATTR_PROTECTED;
}

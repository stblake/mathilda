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

#include <stdio.h>
#include <stdlib.h>

Expr* builtin_algebraicnumberdenominator(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;

    const Expr* a = res->data.function.args[0];
    Expr* out = NULL;
    int r = flint_qqbar_algebraic_number_denominator(a, &out);
    if (r > 0) return out;                   /* success: the denominator */
    if (r < 0) return NULL;                  /* undecided (FLINT compiled out) */

    /* r == 0: not a constant algebraic number — report and stay unevaluated. */
    char* s = expr_to_string((Expr*)a);
    fprintf(stderr,
            "AlgebraicNumberDenominator::nalg: %s is not an explicit "
            "algebraic number.\n", s ? s : "?");
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

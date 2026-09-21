/* algebraicintegerq.c — AlgebraicIntegerQ[x].
 *
 * See algebraicintegerq.h. The exact integrality test (monic primitive minimal
 * polynomial) is delegated to flint_qqbar_algebraic_integer_q; this file maps its
 * tri-state result to True / False / unevaluated.
 */

#include "algebraicintegerq.h"

#include "flint_qqbar.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"

Expr* builtin_algebraicintegerq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;

    int q = flint_qqbar_algebraic_integer_q(res->data.function.args[0]);
    if (q < 0) return NULL;                 /* undecided (FLINT compiled out) */
    return expr_new_symbol(q ? SYM_True : SYM_False);
}

void algebraicintegerq_init(void) {
    symtab_add_builtin("AlgebraicIntegerQ", builtin_algebraicintegerq);
    /* Protected, not Listable: a predicate on the whole expression (a list is
     * not an algebraic integer). */
    symtab_get_def("AlgebraicIntegerQ")->attributes |= ATTR_PROTECTED;
}

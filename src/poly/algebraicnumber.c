/* algebraicnumber.c — AlgebraicNumber[theta, {c0..cn}] canonicalisation.
 *
 * See algebraicnumber.h for the contract. All algebraic-number arithmetic is
 * delegated to flint_qqbar.c (FLINT's exact qqbar engine); this builtin only
 * validates arguments, canonicalises, and guards against re-evaluation churn.
 */

#include "algebraicnumber.h"

#include "flint_qqbar.h"
#include "symtab.h"
#include "attr.h"

Expr* builtin_algebraicnumber(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    const Expr* gen    = res->data.function.args[0];
    const Expr* coeffs = res->data.function.args[1];

    /* Canonical form: reduced AlgebraicNumber[g, {..}] or a rational number.
     * flint_qqbar_algebraic_number declines (NULL) for a non-algebraic
     * generator, malformed coefficients, a degree-cap overflow, or FLINT off. */
    Expr* cand = flint_qqbar_algebraic_number(gen, coeffs);
    if (!cand) return NULL;

    /* Fixpoint guard: if the input is already canonical, leave it untouched so
     * the evaluator does not churn (flint_qqbar_algebraic_number is idempotent,
     * so canon(canon(x)) == canon(x)). */
    if (expr_eq(cand, res)) { expr_free(cand); return NULL; }
    return cand;
}

void algebraicnumber_init(void) {
    symtab_add_builtin("AlgebraicNumber", builtin_algebraicnumber);
    /* NHoldAll: the generator and coefficients must not be numericalised in
     * place (N must reach the dedicated AlgebraicNumber branch, not thread into
     * the args). Numeric behaviour comes from the is_numeric_quantity clause and
     * the numericalize branch, NOT from NumericFunction (whose is_numeric_quantity
     * recursion would fail on the coefficient List). */
    symtab_get_def("AlgebraicNumber")->attributes |= ATTR_NHOLDALL | ATTR_PROTECTED;
}

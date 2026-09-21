/* algebraicnumberpolynomial.c — AlgebraicNumberPolynomial[a, x].
 *
 * See algebraicnumberpolynomial.h for the contract. The coefficient vector of
 * an AlgebraicNumber object is stored in the object itself, so producing the
 * defining polynomial c0 + c1 x + ... + cn x^n is a purely structural read +
 * build: no field arithmetic and no FLINT. The returned Plus/Times/Power tree
 * is canonicalised by the evaluator on its next fixed-point step (zero terms
 * drop, Times[1, x] folds to x, monomials sort by degree).
 */

#include "algebraicnumberpolynomial.h"

#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>

/* True when e is f[...] with head symbol interned as `name` (identity compare). */
static int anp_head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == name;
}

/* An AlgebraicNumber coefficient is an exact rational: Integer, BigInt, or
 * Rational[p, q] (the forms canonicalisation produces and stores). */
static int anp_is_rational_coeff(const Expr* c) {
    return c && (c->type == EXPR_INTEGER || c->type == EXPR_BIGINT ||
                 anp_head_is(c, SYM_Rational));
}

Expr* builtin_algebraicnumberpolynomial(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    const Expr* a = res->data.function.args[0];
    const Expr* x = res->data.function.args[1];

    /* A rational number is the constant polynomial: return it unchanged. */
    if (a->type == EXPR_INTEGER || a->type == EXPR_BIGINT ||
        anp_head_is(a, SYM_Rational))
        return expr_copy((Expr*)a);

    /* AlgebraicNumber[theta, {c0, ..., cn}] with an exact-rational coeff list:
     * build c0 + c1 x + ... + cn x^n. The generator theta is irrelevant here. */
    if (anp_head_is(a, SYM_AlgebraicNumber) && a->data.function.arg_count == 2 &&
        anp_head_is(a->data.function.args[1], SYM_List)) {
        const Expr* coeffs = a->data.function.args[1];
        size_t n = coeffs->data.function.arg_count;

        int all_rational = 1;
        for (size_t i = 0; i < n; i++)
            if (!anp_is_rational_coeff(coeffs->data.function.args[i])) {
                all_rational = 0;
                break;
            }

        if (all_rational) {
            if (n == 0) return expr_new_integer(0);  /* empty list -> 0 */

            Expr** terms = malloc(n * sizeof(Expr*));
            if (!terms) return NULL;
            for (size_t i = 0; i < n; i++) {
                Expr* c = expr_copy(coeffs->data.function.args[i]);
                if (i == 0) {
                    terms[i] = c;                      /* constant term c0 */
                } else if (i == 1) {
                    Expr* t[2] = { c, expr_copy((Expr*)x) };
                    terms[i] = expr_new_function(expr_new_symbol(SYM_Times), t, 2);
                } else {
                    Expr* p[2] = { expr_copy((Expr*)x), expr_new_integer((int64_t)i) };
                    Expr* pw = expr_new_function(expr_new_symbol(SYM_Power), p, 2);
                    Expr* t[2] = { c, pw };
                    terms[i] = expr_new_function(expr_new_symbol(SYM_Times), t, 2);
                }
            }
            Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, n);
            free(terms);
            return sum;
        }
    }

    /* Not a valid AlgebraicNumber object or rational: warn and decline. */
    char* s = expr_to_string((Expr*)a);
    fprintf(stderr,
            "AlgebraicNumberPolynomial::naobj: %s is not a valid "
            "AlgebraicNumber object.\n", s ? s : "?");
    free(s);
    return NULL;
}

void algebraicnumberpolynomial_init(void) {
    symtab_add_builtin("AlgebraicNumberPolynomial",
                       builtin_algebraicnumberpolynomial);
    symtab_get_def("AlgebraicNumberPolynomial")->attributes |=
        ATTR_LISTABLE | ATTR_PROTECTED;
}

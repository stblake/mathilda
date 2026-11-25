/* definite_minors.c - definiteness of a matrix of exact numeric constants.
 *
 * PositiveDefiniteMatrixQ / NegativeDefiniteMatrixQ decide an all-rational
 * matrix exactly in FLINT and an approximate one by Cholesky in doubles. An
 * entry such as Pi or Sqrt[2] fits neither: rounding it to a double would make
 * the answer depend on rounding, and FLINT has no such field. Sylvester's
 * criterion does: the Hermitian part H = (m + m^H)/2 is positive definite iff
 * every leading principal minor is positive, and the sign of an exact closed
 * form is certified by Positive (a rigorous ball evaluation, else an exact
 * zero test -- never a rounded value). An undecided sign leaves the predicate
 * to answer False, the convention for "not provably definite".
 *
 * Memory (SPEC section 4): `m` is borrowed; nothing is returned but an int.
 */

#include "linalg.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include <stdbool.h>
#include <stdlib.h>

/* Leading minors are symbolic determinants: past this size their cost is not
 * worth paying for a predicate. */
#define DEFINITE_MINORS_MAX_N 8

static bool has_inexact(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_inexact(e->data.function.args[i])) return true;
    return false;
}

bool linalg_is_exact_constant(const Expr* e) {
    return e && e->type != EXPR_INTEGER && e->type != EXPR_BIGINT &&
           !has_inexact(e) && expr_is_numeric_quantity(e);
}

static Expr* call1(const char* h, Expr* a) {
    return eval_and_free(expr_new_function(expr_new_symbol(h), &a, 1));
}

int linalg_exact_definite_minors(Expr* m, int n, int negative) {
    if (n < 1 || n > DEFINITE_MINORS_MAX_N) return -1;
    /* H = (m + ConjugateTranspose[m]) / 2 */
    Expr* ct = call1(SYM_ConjugateTranspose, expr_copy(m));
    Expr* sum = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(m), ct }, 2));
    Expr* half = eval_and_free(expr_new_function(expr_new_symbol(SYM_Rational),
                    (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2));
    Expr* H = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ half, sum }, 2));
    if (H->type != EXPR_FUNCTION || (int)H->data.function.arg_count != n) {
        expr_free(H);
        return -1;
    }
    int verdict = 1;
    for (int k = 1; k <= n && verdict == 1; k++) {
        /* the k x k leading block */
        Expr** rows = malloc(sizeof(Expr*) * (size_t)k);
        for (int i = 0; i < k; i++) {
            Expr** el = malloc(sizeof(Expr*) * (size_t)k);
            const Expr* hr = H->data.function.args[i];
            for (int j = 0; j < k; j++) el[j] = expr_copy(hr->data.function.args[j]);
            rows[i] = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)k);
            free(el);
        }
        Expr* blk = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)k);
        free(rows);
        Expr* d = call1(SYM_Det, blk);
        if (negative && (k % 2 == 1))
            d = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), d }, 2));
        Expr* s = call1(SYM_Positive, d);
        if (s->type == EXPR_SYMBOL && s->data.symbol.name == SYM_True) verdict = 1;
        else if (s->type == EXPR_SYMBOL && s->data.symbol.name == SYM_False) verdict = 0;
        else verdict = -1;
        expr_free(s);
    }
    expr_free(H);
    return verdict;
}

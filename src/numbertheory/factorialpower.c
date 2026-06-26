/* factorialpower.c -- FactorialPower[].
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

/* FactorialPower[n, k] -- the falling factorial (Pochhammer-style).
 *
 *   FactorialPower[n, k] = n * (n - 1) * (n - 2) * ... * (n - k + 1)
 *                        = n! / (n - k)!     (for non-negative integer k)
 *                        = 1                  (k == 0)
 *
 * Concrete numeric paths:
 *   - Both arguments integer with k >= 0: compute the product (GMP for
 *     k > 30 to avoid int64 overflow).
 *   - Integer n with negative integer k: undefined / no closed-form
 *     elementary identity in MMA's symbolic core; leave unevaluated.
 *
 * Symbolic n with concrete non-negative integer k expands to a product
 * of (n - i) factors so downstream Expand / D can act on it -- e.g.
 *   FactorialPower[n, 3] -> n (n - 1) (n - 2)
 * Mathematically equal under Mathematica's definition and lets the
 * surrounding evaluator collapse it via Expand[] when the user wants.
 *
 * Anything else (k symbolic, or both arguments symbolic): NULL --
 * stays as `FactorialPower[n, k]` until more information is supplied.
 */
Expr* builtin_factorialpower(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* n_arg = res->data.function.args[0];
    Expr* k_arg = res->data.function.args[1];

    /* k must be a concrete non-negative integer for any reduction. */
    if (k_arg->type != EXPR_INTEGER) return NULL;
    int64_t k = k_arg->data.integer;
    if (k < 0) return NULL;
    if (k == 0) return expr_new_integer(1);

    /* Both n and k integer: compute exact product in GMP. */
    if (n_arg->type == EXPR_INTEGER || n_arg->type == EXPR_BIGINT) {
        mpz_t n_z, term, product;
        mpz_init(n_z);
        if (n_arg->type == EXPR_INTEGER) mpz_set_si(n_z, n_arg->data.integer);
        else                              mpz_set(n_z, n_arg->data.bigint);
        mpz_init_set_ui(product, 1);
        mpz_init(term);
        for (int64_t i = 0; i < k; i++) {
            mpz_set(term, n_z);
            mpz_sub_ui(term, term, (unsigned long)i);
            mpz_mul(product, product, term);
        }
        mpz_clear(n_z);
        mpz_clear(term);
        Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(product));
        mpz_clear(product);
        return result;
    }

    /* Symbolic n, k a small concrete integer (<= 32 to keep tree size
     * bounded). Expand to Times[n, n-1, ..., n-k+1]. */
    if (k <= 32) {
        Expr** factors = malloc(sizeof(Expr*) * (size_t)k);
        for (int64_t i = 0; i < k; i++) {
            if (i == 0) {
                factors[i] = expr_copy(n_arg);
            } else {
                Expr* shift = expr_new_integer(-i);
                factors[i] = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(n_arg), shift }, 2);
            }
        }
        Expr* product = expr_new_function(expr_new_symbol(SYM_Times), factors, (size_t)k);
        free(factors);
        return eval_and_free(product);
    }
    return NULL;
}

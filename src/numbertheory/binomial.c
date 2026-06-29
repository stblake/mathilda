/* binomial.c -- Binomial[] and binomial_polynomial helper.
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

/* Binomial[n, k] expanded as a falling-factorial polynomial for a small
 * non-negative integer k and arbitrary n (symbolic, rational, complex,
 * real, integer, ...).  Yields
 *
 *     Binomial[n, k] = n (n - 1) (n - 2) ... (n - k + 1) / k!
 *
 * Returns a fully-evaluated expression.  Used for two cases:
 *   - concrete non-negative integer k <= 32, n anything non-integer:
 *     Binomial[n, 4]  -> n (n-1) (n-2) (n-3) / 24
 *     Binomial[1+I,5] -> -1/12 - I/12 (Times/Plus fold complex factors)
 *   - the symmetric-reduction case where Subtract[n, m] collapses to a
 *     small non-negative integer k: Binomial[9/2, 7/2] -> Binomial[9/2,1]
 *     -> 9/2; Binomial[n, n-1] -> Binomial[n, 1] -> n.
 *
 * Cap of 32 mirrors FactorialPower so a careless Binomial[n, 1000] does
 * not blow up the tree size.  Returns NULL for k outside [0, 32]. */
static Expr* binomial_polynomial(Expr* n, int64_t k) {
    if (k < 0 || k > 32) return NULL;
    if (k == 0) return expr_new_integer(1);
    if (k == 1) return expr_copy(n);

    /* numerator factors: n, n-1, ..., n-k+1, plus a trailing 1/k! coefficient. */
    Expr** factors = malloc(sizeof(Expr*) * (size_t)(k + 1));
    for (int64_t i = 0; i < k; i++) {
        if (i == 0) {
            factors[i] = expr_copy(n);
        } else {
            Expr* shift = expr_new_integer(-i);
            factors[i] = expr_new_function(expr_new_symbol(SYM_Plus),
                (Expr*[]){ expr_copy(n), shift }, 2);
        }
    }

    /* 1 / k! coefficient.  k! is computed in GMP so values for k > 20
     * (where k! overflows int64) still round-trip exactly through a
     * Rational[1, BigInt] which the evaluator's is_rational_like folders
     * recognise. */
    mpz_t fact_z;
    mpz_init(fact_z);
    mpz_fac_ui(fact_z, (unsigned long)k);
    Expr* fact_expr = expr_bigint_normalize(expr_new_bigint_from_mpz(fact_z));
    mpz_clear(fact_z);

    Expr* inv;
    if (fact_expr->type == EXPR_INTEGER) {
        inv = make_rational(1, fact_expr->data.integer);
        expr_free(fact_expr);
    } else {
        Expr* one = expr_new_integer(1);
        inv = expr_new_function(expr_new_symbol(SYM_Rational),
            (Expr*[]){ one, fact_expr }, 2);
    }
    factors[k] = inv;

    Expr* product = expr_new_function(expr_new_symbol(SYM_Times),
        factors, (size_t)(k + 1));
    free(factors);
    return eval_and_free(product);
}

/* Binomial[n, m] -- generalised binomial coefficient
 *
 *     Binomial[n, m] = Gamma[n+1] / (Gamma[m+1] * Gamma[n - m + 1])
 *
 * Evaluation strategy (first matching rule wins):
 *   1. Both args concrete integer-like: exact GMP via mpz_bin_ui, with
 *      Mathematica's Pascal extension Binomial[n, m] = (-1)^m Binomial[
 *      m - n - 1, m] for n < 0, m >= 0.
 *   2. Either arg machine-precision real: tgamma path -> EXPR_REAL.
 *   3. Symmetric identity: if Subtract[n, m] simplifies to a non-negative
 *      integer k <= 32, expand Binomial[n, k] as a falling-factorial
 *      polynomial.  This catches Binomial[n, n-1] -> n,
 *      Binomial[9/2, 7/2] -> 9/2, Binomial[n+1, n-1] -> n (n+1) / 2, ...
 *   4. Concrete non-negative integer m <= 32 with anything-n: expand
 *      Binomial[n, m] as a falling-factorial polynomial.  Handles
 *      Binomial[n, 4] -> n(n-1)(n-2)(n-3)/24 and lets the Times/Plus
 *      folders simplify complex / rational n (Binomial[1+I, 5]
 *      -> -1/12 - I/12). */
Expr* builtin_binomial(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* arg_n = res->data.function.args[0];
    Expr* arg_m = res->data.function.args[1];

    /* --- (1) exact integer / integer path. --- */
    if (expr_is_integer_like(arg_n) && expr_is_integer_like(arg_m)) {
        mpz_t n, m;
        expr_to_mpz(arg_n, n);
        expr_to_mpz(arg_m, m);

        if (mpz_sgn(m) < 0) {
            mpz_clears(n, m, NULL);
            return expr_new_integer(0);
        }
        if (mpz_sgn(n) >= 0 && mpz_cmp(m, n) > 0) {
            mpz_clears(n, m, NULL);
            return expr_new_integer(0);
        }
        /* mpz_bin_ui needs m to fit in unsigned long; pathologically
         * large m has no finite expansion we'd want anyway. */
        if (!mpz_fits_ulong_p(m)) {
            mpz_clears(n, m, NULL);
            return NULL;
        }
        unsigned long mk = mpz_get_ui(m);
        mpz_t r;
        mpz_init(r);
        if (mpz_sgn(n) < 0) {
            /* Pascal extension for negative n. */
            mpz_t shifted;
            mpz_init(shifted);
            mpz_neg(shifted, n);
            mpz_add(shifted, shifted, m);
            mpz_sub_ui(shifted, shifted, 1);
            mpz_bin_ui(r, shifted, mk);
            if (mk & 1UL) mpz_neg(r, r);
            mpz_clear(shifted);
        } else {
            mpz_bin_ui(r, n, mk);
        }
        Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
        mpz_clears(n, m, r, NULL);
        return out;
    }

    /* --- (2) machine real path via tgamma. --- */
    if (arg_n->type == EXPR_REAL || arg_m->type == EXPR_REAL) {
        double nv = (arg_n->type == EXPR_REAL)    ? arg_n->data.real :
                    (arg_n->type == EXPR_INTEGER) ? (double)arg_n->data.integer :
                    (arg_n->type == EXPR_BIGINT)  ? mpz_get_d(arg_n->data.bigint) : 0.0;
        double mv = (arg_m->type == EXPR_REAL)    ? arg_m->data.real :
                    (arg_m->type == EXPR_INTEGER) ? (double)arg_m->data.integer :
                    (arg_m->type == EXPR_BIGINT)  ? mpz_get_d(arg_m->data.bigint) : 0.0;
        return expr_new_real(tgamma(nv + 1.0) / (tgamma(mv + 1.0) * tgamma(nv - mv + 1.0)));
    }

    /* --- (3) symmetric identity: reduce when n - m is a small non-neg int.
     *
     * Warnings are muted around the Subtract evaluation: when n / m happen
     * to involve a 1/0 partial that lives only inside this exploratory
     * difference, the user is not subtracting -- we are -- so the
     * Power::infy diagnostic would be spurious. --- */
    {
        arith_warnings_mute_push();
        Expr* diff = expr_new_function(expr_new_symbol(SYM_Subtract),
            (Expr*[]){ expr_copy(arg_n), expr_copy(arg_m) }, 2);
        Expr* diff_eval = eval_and_free(diff);
        arith_warnings_mute_pop();

        int64_t k = -1;
        if (diff_eval && diff_eval->type == EXPR_INTEGER &&
            diff_eval->data.integer >= 0 && diff_eval->data.integer <= 32) {
            k = diff_eval->data.integer;
        }
        expr_free(diff_eval);

        if (k >= 0) {
            Expr* poly = binomial_polynomial(arg_n, k);
            if (poly) return poly;
        }
    }

    /* --- (4) concrete non-negative integer m: falling-factorial expand. --- */
    if (arg_m->type == EXPR_INTEGER && arg_m->data.integer >= 0 &&
        arg_m->data.integer <= 32) {
        return binomial_polynomial(arg_n, arg_m->data.integer);
    }

    return NULL;
}

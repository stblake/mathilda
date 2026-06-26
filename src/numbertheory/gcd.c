/* gcd.c -- GCD[].
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

/* GCD over rational-like args folded entirely in GMP: gcd(a/b, c/d) =
 * gcd(a, c) / lcm(b, d). Used both for the BigInt-bearing case and as the
 * int64 fast-path overflow fallback. Every arg must be rational-like. */
static Expr* gcd_rational_mpz_fold(Expr* res, size_t count) {
    mpz_t running_n, running_d, n, d;
    mpz_init_set_ui(running_n, 0);
    mpz_init_set_ui(running_d, 1);
    mpz_inits(n, d, NULL);
    for (size_t i = 0; i < count; i++) {
        rational_like_to_mpz_pair(res->data.function.args[i], n, d);
        mpz_abs(n, n);
        mpz_abs(d, d);
        if (i == 0) {
            mpz_set(running_n, n);
            mpz_set(running_d, d);
        } else {
            mpz_gcd(running_n, running_n, n);
            mpz_lcm(running_d, running_d, d);
        }
    }
    Expr* result = mpz_pair_to_rational_expr(running_n, running_d);
    mpz_clears(running_n, running_d, n, d, NULL);
    return result;
}
Expr* builtin_gcd(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t count = res->data.function.arg_count;
    if (count == 0) return expr_new_integer(0);
    if (count == 1) return single_arg_abs_or_copy(res->data.function.args[0]);

    // Single pass: detect any bigint while confirming all args are integer-like
    bool any_bigint = false, all_integer_like = true;
    bool all_rational_like = true, any_rational_needs_bigint = false;
    for (size_t i = 0; i < count; i++) {
        Expr* arg = res->data.function.args[i];
        if (!expr_is_integer_like(arg)) all_integer_like = false;
        else if (arg->type == EXPR_BIGINT) any_bigint = true;
        if (!is_rational_like(arg)) all_rational_like = false;
        else if (rational_like_needs_bigint(arg)) any_rational_needs_bigint = true;
    }

    if (all_integer_like && any_bigint) {
        mpz_t running, tmp;
        mpz_init_set_ui(running, 0);
        for (size_t i = 0; i < count; i++) {
            expr_to_mpz(res->data.function.args[i], tmp);
            mpz_abs(tmp, tmp);
            mpz_gcd(running, running, tmp);
            mpz_clear(tmp);
        }
        Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(running));
        mpz_clear(running);
        return result;
    }

    /* Rational-bigint fold. Triggers when every arg is rational-like and
     * at least one carries a BigInt component (bare BigInt or BigInt
     * inside Rational). gcd(a/b, c/d) = gcd(a, c) / lcm(b, d) folded
     * pairwise. */
    if (all_rational_like && any_rational_needs_bigint) {
        return gcd_rational_mpz_fold(res, count);
    }

    int64_t running_n = 0;
    int64_t running_d = 1;
    bool overflow = false;

    for (size_t i = 0; i < count; i++) {
        int64_t n, d;
        if (!is_rational(res->data.function.args[i], &n, &d)) {
            return NULL;
        }
        if (i == 0) {
            running_n = llabs(n);
            running_d = llabs(d);
        } else {
            running_n = gcd(running_n, n);
            /* gcd(a/b, c/d) denominator is lcm(b, d), which can overflow. */
            running_d = lcm_checked(running_d, d, &overflow);
            if (overflow) return gcd_rational_mpz_fold(res, count);
        }
    }

    return make_rational(running_n, running_d);
}

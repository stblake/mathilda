/* lcm.c -- LCM[].
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

/* LCM over rational-like args folded entirely in GMP: lcm(a/b, c/d) =
 * lcm(a, c) / gcd(b, d), with a zero numerator zeroing the running LCM.
 * Used for the BigInt-bearing case and the int64 overflow fallback. */
static Expr* lcm_rational_mpz_fold(Expr* res, size_t count) {
    mpz_t running_n, running_d, n, d;
    mpz_init_set_ui(running_n, 0);
    mpz_init_set_ui(running_d, 0);
    mpz_inits(n, d, NULL);
    for (size_t i = 0; i < count; i++) {
        rational_like_to_mpz_pair(res->data.function.args[i], n, d);
        mpz_abs(n, n);
        mpz_abs(d, d);
        if (i == 0) {
            mpz_set(running_n, n);
            mpz_set(running_d, d);
        } else if (mpz_sgn(running_n) == 0 || mpz_sgn(n) == 0) {
            mpz_set_ui(running_n, 0);
            mpz_set_ui(running_d, 1);
        } else {
            mpz_lcm(running_n, running_n, n);
            mpz_gcd(running_d, running_d, d);
        }
    }
    Expr* result = mpz_pair_to_rational_expr(running_n, running_d);
    mpz_clears(running_n, running_d, n, d, NULL);
    return result;
}
Expr* builtin_lcm(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t count = res->data.function.arg_count;
    if (count == 0) return expr_new_integer(1);
    if (count == 1) return single_arg_abs_or_copy(res->data.function.args[0]);

    /* Mirror builtin_gcd: when every arg is integer-like and any one is a
     * bigint, fold over GMP so results past int64 don't overflow. */
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
        mpz_init_set_ui(running, 1);
        for (size_t i = 0; i < count; i++) {
            expr_to_mpz(res->data.function.args[i], tmp);
            mpz_abs(tmp, tmp);
            mpz_lcm(running, running, tmp);
            mpz_clear(tmp);
            /* mpz_lcm(_, _, 0) yields 0; any further LCMs stay 0. */
            if (mpz_sgn(running) == 0) break;
        }
        Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(running));
        mpz_clear(running);
        return result;
    }

    /* Rational-bigint fold. lcm(a/b, c/d) = lcm(a, c) / gcd(b, d) folded
     * pairwise. A zero numerator in any term zeroes the running lcm. */
    if (all_rational_like && any_rational_needs_bigint) {
        return lcm_rational_mpz_fold(res, count);
    }

    int64_t running_n = 0;
    int64_t running_d = 0;
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
            if (running_n == 0 || n == 0) {
                running_n = 0;
                running_d = 1;
            } else {
                running_n = lcm_checked(running_n, n, &overflow);
                running_d = gcd(running_d, d);
                if (overflow) return lcm_rational_mpz_fold(res, count);
            }
        }
    }

    return make_rational(running_n, running_d);
}

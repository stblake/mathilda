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
#include "message.h"   /* mth_message_gated: Quiet/Check funnel */
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

/* LCM over the Gaussian rationals, used when some argument is a non-real
 * exact complex number: lcm(w_1, w_2, ...) / gcd(d_1, d_2, ...) for arguments
 * w_i/d_i, normalised to the first-quadrant associate -- LCM[3 + I, 2] is
 * 2 + 4 I and LCM[2 + I, 2 - I] is 5. Any zero argument
 * gives 0. NULL if some argument is not a Gaussian rational. */
static Expr* lcm_gaussian(Expr* res, size_t count) {
    mpz_t wr, wi, d, xr, xi, xd, lr, li;
    mpz_inits(wr, wi, d, xr, xi, xd, lr, li, NULL);
    bool ok = true;
    for (size_t i = 0; i < count && ok; i++) {
        if (!gq_from_expr(res->data.function.args[i], xr, xi, xd)) { ok = false; break; }
        if (i == 0) {
            mpz_set(wr, xr); mpz_set(wi, xi); mpz_set(d, xd);
            gi_normalize(wr, wi, NULL);
        } else {
            gi_lcm(lr, li, wr, wi, xr, xi);
            mpz_swap(wr, lr); mpz_swap(wi, li);
            mpz_gcd(d, d, xd);
        }
    }
    if (ok && mpz_sgn(wr) == 0 && mpz_sgn(wi) == 0) mpz_set_ui(d, 1);
    Expr* out = ok ? gq_to_expr(wr, wi, d) : NULL;
    mpz_clears(wr, wi, d, xr, xi, xd, lr, li, NULL);
    return out;
}

Expr* builtin_lcm(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t count = res->data.function.arg_count;
    /* LCM::argm -- unlike GCD[] (0), LCM[] has no value. */
    if (count == 0) {
        mth_message_gated(g_arith_warnings_muted, "LCM", "argm",
                          "LCM called with 0 arguments; 1 or more arguments are expected.");
        return NULL;
    }

    /* LCM::exact -- an inexact number has no lcm (LCM[2.5, 5]). */
    Expr* inexact = nt_first_inexact_arg(res);
    if (inexact) return nt_emit_exact("LCM", inexact, res);

    /* Any non-real exact complex argument switches to Z[i]. */
    for (size_t i = 0; i < count; i++)
        if (is_gaussian_rational(res->data.function.args[i]))
            return lcm_gaussian(res, count);

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
        /* INT64_MIN has no int64 magnitude (llabs leaves it negative, and
         * gcd(-2^63, 0) = 2^63 does not fit): GCD[-2^63] was -2^63. */
        if (n == INT64_MIN || d == INT64_MIN) return lcm_rational_mpz_fold(res, count);
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

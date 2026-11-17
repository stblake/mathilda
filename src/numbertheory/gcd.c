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
#include "message.h"   /* mth_message_gated: Quiet/Check funnel */
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

/* GCD over the Gaussian rationals, used when some argument is a non-real
 * exact complex number: writing each argument as w/d (w in Z[i], d a positive
 * integer), the gcd is gcd(w_1, w_2, ...) / lcm(d_1, d_2, ...) with the
 * Gaussian gcd normalised to its first-quadrant associate -- GCD[3 + I, 2] is
 * 1 + I and GCD[1/2 + I/3, 2] is 1/6. NULL if some argument
 * is not a Gaussian rational (GCD[2, 3 + I, x] stays unevaluated). */
static Expr* gcd_gaussian(Expr* res, size_t count) {
    mpz_t wr, wi, d, xr, xi, xd, gr, gi;
    mpz_inits(wr, wi, d, xr, xi, xd, gr, gi, NULL);
    mpz_set_ui(d, 1);
    bool ok = true;
    for (size_t i = 0; i < count && ok; i++) {
        if (!gq_from_expr(res->data.function.args[i], xr, xi, xd)) { ok = false; break; }
        gi_gcd(gr, gi, wr, wi, xr, xi);
        mpz_swap(wr, gr); mpz_swap(wi, gi);
        mpz_lcm(d, d, xd);
    }
    Expr* out = ok ? gq_to_expr(wr, wi, d) : NULL;
    mpz_clears(wr, wi, d, xr, xi, xd, gr, gi, NULL);
    return out;
}

Expr* builtin_gcd(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t count = res->data.function.arg_count;
    if (count == 0) return expr_new_integer(0);

    /* GCD::exact -- an inexact number has no gcd (GCD[2.5, 5]). */
    Expr* inexact = nt_first_inexact_arg(res);
    if (inexact) return nt_emit_exact("GCD", inexact, res);

    /* Any non-real exact complex argument switches to Z[i] (a lone one is
     * normalised: GCD[-1 - I] is 1 + I, GCD[2 I] is 2). */
    for (size_t i = 0; i < count; i++)
        if (is_gaussian_rational(res->data.function.args[i]))
            return gcd_gaussian(res, count);

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
        /* INT64_MIN has no int64 magnitude (llabs leaves it negative, and
         * gcd(-2^63, 0) = 2^63 does not fit): GCD[-2^63] was -2^63. */
        if (n == INT64_MIN || d == INT64_MIN) return gcd_rational_mpz_fold(res, count);
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

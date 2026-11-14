/* extendedgcd.c -- ExtendedGCD[].
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
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

/* ----------------------------------------------------------------------
 * ExtendedGCD[n1, n2, ...] -> {g, {r1, r2, ...}} where g == GCD[n1, ...]
 * and g == r1 n1 + r2 n2 + ... (a Bezout / multi-argument extended GCD).
 *
 * Integer-only: machine ints and GMP bigints are both accepted (the
 * int64 fast path auto-promotes through expr_to_mpz, and results demote
 * back to EXPR_INTEGER via expr_bigint_normalize when they fit).
 *
 * The cofactors are produced by folding GMP's mpz_gcdext pairwise:
 *
 *     gcd(running_g, a_i) = s * running_g + t * a_i
 *
 * Each fold scales every coefficient accumulated so far by s and appends
 * t for the new argument. running_g starts at 0 (so the first step yields
 * |a_0| with cofactor sign(a_0)) and stays non-negative throughout — GMP
 * normalises the gcd non-negative — so the leading g is exactly GCD and
 * the Bezout identity holds. Matches Mathematica's sign convention. */

/* ExtendedGCD::exact -- an inexact (Real / MPFR) argument was supplied. */
static Expr* egcd_emit_exact(Expr* arg, Expr* res) {
    if (!g_arith_warnings_muted) {
        char* arg_str = expr_to_string(arg);
        char* call_str = expr_to_string(res);
        mth_message("ExtendedGCD", "exact", "Argument %s in %s is not an exact number.", arg_str ? arg_str : "?", call_str ? call_str : "?");
        free(arg_str);
        free(call_str);
    }
    return NULL;
}

/* ExtendedGCD::egcd -- every argument is exact but at least one is a
 * non-integer rational. */
static Expr* egcd_emit_egcd(Expr* res) {
    if (!g_arith_warnings_muted) {
        char* call_str = expr_to_string(res);
        mth_message("ExtendedGCD", "egcd", "Arguments in %s should be integers.", call_str ? call_str : "?");
        free(call_str);
    }
    return NULL;
}

/* {g, {c1, c2, ...}} from mpz values. */
static Expr* egcd_result(Expr* g, Expr** coeffs, size_t count) {
    Expr* coeff_list = expr_new_function(expr_new_symbol(SYM_List), coeffs, count);
    Expr* pair[2] = { g, coeff_list };
    return expr_new_function(expr_new_symbol(SYM_List), pair, 2);
}

static Expr* egcd_gaussian_value(const mpz_t re, const mpz_t im) {
    mpz_t one;
    mpz_init_set_ui(one, 1);
    Expr* e = gq_to_expr(re, im, one);
    mpz_clear(one);
    return e;
}

/* ExtendedGCD over Z[i], used when some argument is a non-real Gaussian
 * integer. Folds gi_xgcd from the first argument (normalised, with the unit
 * that normalised it as its coefficient), scaling the accumulated coefficients
 * by each step's s -- the integer fold's shape. Gives the canonical
 * coefficients, e.g. ExtendedGCD[3 + I, 2] is {1 + I, {-I, 2 I}}; see gi_xgcd
 * for where the two can differ. Every argument must be a Gaussian integer. */
static Expr* egcd_gaussian(Expr* res, size_t count) {
    mpz_t gr, gi, sr, si, tr, ti, ar, ai, nr, ni, pr, pim;
    mpz_inits(gr, gi, sr, si, tr, ti, ar, ai, nr, ni, pr, pim, NULL);
    mpz_t* cr = malloc(count * sizeof(mpz_t));
    mpz_t* ci = malloc(count * sizeof(mpz_t));
    for (size_t i = 0; i < count; i++) { mpz_init(cr[i]); mpz_init(ci[i]); }
    for (size_t i = 0; i < count; i++) {
        mpz_t ta, tb;
        if (!df_to_gaussian(res->data.function.args[i], ta, tb)) continue;  /* checked by caller */
        mpz_set(ar, ta); mpz_set(ai, tb);
        mpz_clears(ta, tb, NULL);
        if (i == 0) {
            mpz_set(gr, ar); mpz_set(gi, ai);
            int k;
            gi_normalize(gr, gi, &k);
            mpz_set_ui(cr[0], 1); mpz_set_ui(ci[0], 0);
            gi_rotate_by(cr[0], ci[0], k);
            continue;
        }
        gi_xgcd(nr, ni, sr, si, tr, ti, gr, gi, ar, ai);
        for (size_t j = 0; j < i; j++) {             /* c_j <- c_j * s */
            mpz_mul(pr, cr[j], sr); mpz_submul(pr, ci[j], si);
            mpz_mul(pim, cr[j], si); mpz_addmul(pim, ci[j], sr);
            mpz_swap(cr[j], pr); mpz_swap(ci[j], pim);
        }
        mpz_set(cr[i], tr); mpz_set(ci[i], ti);
        mpz_set(gr, nr); mpz_set(gi, ni);
    }
    Expr** coeffs = malloc(count * sizeof(Expr*));
    for (size_t i = 0; i < count; i++) {
        coeffs[i] = egcd_gaussian_value(cr[i], ci[i]);
        mpz_clears(cr[i], ci[i], NULL);
    }
    free(cr); free(ci);
    Expr* out = egcd_result(egcd_gaussian_value(gr, gi), coeffs, count);
    free(coeffs);
    mpz_clears(gr, gi, sr, si, tr, ti, ar, ai, nr, ni, pr, pim, NULL);
    return out;
}

/* The two-machine-integer case without GMP. Reproduces mpz_gcdext's unique
 * cofactors -- |s| < |b|/(2g) and |t| < |a|/(2g), with s = 0, t = sgn(b) when
 * |a| == |b|, and s = sgn(a) when b == 0 or |b| == 2g, t = sgn(b) when a == 0
 * or |a| == 2g -- by running the extended Euclidean algorithm on |a|, |b| and
 * then reducing (s, t) into that window. Declines (false) when an operand is
 * INT64_MIN; the GMP fold then answers. */
static bool egcd_i64(int64_t a, int64_t b, int64_t* g_out, int64_t* s_out, int64_t* t_out) {
    if (a == INT64_MIN || b == INT64_MIN) return false;
    /* Fold shape of the GMP path: gcdext(0, a) gives (|a|, 0, sgn a); then
     * gcdext(|a|, b) gives (g, s, t), and the cofactors are (sgn(a) s, t). */
    int64_t A = a < 0 ? -a : a;
    int64_t sa = (a > 0) - (a < 0);
    int64_t B = b < 0 ? -b : b;
    int64_t sb = (b > 0) - (b < 0);
    /* mpz_gcdext(A, b) with A >= 0. */
    int64_t g, s, t;
    if (A == 0 && B == 0) { g = 0; s = 0; t = 0; }
    else if (A == B) { g = A; s = 0; t = sb; }
    else if (B == 0) { g = A; s = (A > 0); t = 0; }
    else if (A == 0) { g = B; s = 0; t = sb; }
    else {
        int64_t r0 = A, r1 = B, s0 = 1, s1 = 0, t0 = 0, t1 = 1;
        while (r1 != 0) {
            int64_t q = r0 / r1, tmp;
            tmp = r0 - q * r1; r0 = r1; r1 = tmp;
            tmp = s0 - q * s1; s0 = s1; s1 = tmp;
            tmp = t0 - q * t1; t0 = t1; t1 = tmp;
        }
        g = r0; s = s0; t = t0 * sb;             /* g = s A + t0 B = s A + t b */
        /* Reduce s into |s| < B/(2g) (and t with it); exceptional |B| == 2g
         * gives s = sgn(A) = 1. */
        int64_t bg = B / g, ag = A / g;
        if (B == 2 * g) { s = 1; t = (g - A) / b; }
        else if (A == 2 * g) { t = sb; s = (g - b * t) / A; }
        else {
            while (2 * s >= bg)  { s -= bg; t += ag * sb; }
            while (2 * s <= -bg) { s += bg; t -= ag * sb; }
        }
    }
    *g_out = g; *s_out = sa * s; *t_out = t;
    return true;
}

Expr* builtin_extendedgcd(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t count = res->data.function.arg_count;

    /* ExtendedGCD[] -> {0, {}} */
    if (count == 0) {
        Expr* g = expr_new_integer(0);
        Expr* empty = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        Expr* pair[2] = { g, empty };
        return expr_new_function(expr_new_symbol(SYM_List), pair, 2);
    }

    /* Classify the arguments in a single pass. */
    bool all_integer = true, all_exact = true, any_inexact = false;
    Expr* first_inexact = NULL;
    for (size_t i = 0; i < count; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type == EXPR_REAL
#ifdef USE_MPFR
            || arg->type == EXPR_MPFR
#endif
            ) {
            all_integer = false;
            all_exact = false;
            if (!any_inexact) { any_inexact = true; first_inexact = arg; }
        } else if (expr_is_integer_like(arg)) {
            /* integer: ok */
        } else if (is_rational_like(arg)) {
            all_integer = false;   /* exact, but not an integer */
        } else {
            all_integer = false;   /* symbolic / non-numeric */
            all_exact = false;
        }
    }

    if (any_inexact) return egcd_emit_exact(first_inexact, res);

    /* Z[i]: any non-real Gaussian argument. Non-integer Gaussian rationals
     * draw ::egcd, as their real counterparts do. */
    bool any_gauss = false, all_gauss_int = true;
    for (size_t i = 0; i < count; i++) {
        Expr* arg = res->data.function.args[i];
        if (is_gaussian_rational(arg)) {
            any_gauss = true;
            if (!is_gaussian_integer(arg)) all_gauss_int = false;
        } else if (!expr_is_integer_like(arg)) {
            all_gauss_int = false;
        }
    }
    if (any_gauss) {
        if (all_gauss_int) return egcd_gaussian(res, count);
        bool all_exact_num = true;
        for (size_t i = 0; i < count; i++) {
            Expr* arg = res->data.function.args[i];
            if (!is_gaussian_rational(arg) && !is_rational_like(arg)) all_exact_num = false;
        }
        return all_exact_num ? egcd_emit_egcd(res) : NULL;
    }
    if (!all_integer) {
        /* All-exact-but-non-integer -> diagnose; anything symbolic stays
         * unevaluated silently so rules / symbolic flow still apply. */
        if (all_exact) return egcd_emit_egcd(res);
        return NULL;
    }

    /* The common call, two machine integers, needs no GMP. */
    if (count == 2 && res->data.function.args[0]->type == EXPR_INTEGER &&
        res->data.function.args[1]->type == EXPR_INTEGER) {
        int64_t g, s1, t1;
        if (egcd_i64(res->data.function.args[0]->data.integer,
                     res->data.function.args[1]->data.integer, &g, &s1, &t1)) {
            Expr* coeffs[2] = { expr_new_integer(s1), expr_new_integer(t1) };
            return egcd_result(expr_new_integer(g), coeffs, 2);
        }
    }

    /* --- Integer fold via the extended Euclidean algorithm (GMP). --- */
    mpz_t running_g, new_g, s, t;
    mpz_inits(running_g, new_g, s, t, NULL);
    mpz_set_ui(running_g, 0);

    /* coeffs[i] accumulates the cofactor of argument i. */
    mpz_t* coeffs = malloc(count * sizeof(mpz_t));
    if (!coeffs) {
        mpz_clears(running_g, new_g, s, t, NULL);
        return NULL;
    }

    for (size_t i = 0; i < count; i++) {
        mpz_t a;
        expr_to_mpz(res->data.function.args[i], a); /* inits a */
        mpz_gcdext(new_g, s, t, running_g, a);
        /* Scale previously accumulated cofactors by s, then append t. */
        for (size_t j = 0; j < i; j++) {
            mpz_mul(coeffs[j], coeffs[j], s);
        }
        mpz_init_set(coeffs[i], t);
        mpz_set(running_g, new_g);
        mpz_clear(a);
    }

    Expr* g_expr = expr_bigint_normalize(expr_new_bigint_from_mpz(running_g));
    Expr** coeff_args = malloc(count * sizeof(Expr*));
    if (!coeff_args) {
        /* Out of memory: clean up and bail (leave call unevaluated). */
        expr_free(g_expr);
        for (size_t i = 0; i < count; i++) mpz_clear(coeffs[i]);
        free(coeffs);
        mpz_clears(running_g, new_g, s, t, NULL);
        return NULL;
    }
    for (size_t i = 0; i < count; i++) {
        coeff_args[i] = expr_bigint_normalize(expr_new_bigint_from_mpz(coeffs[i]));
        mpz_clear(coeffs[i]);
    }
    free(coeffs);
    mpz_clears(running_g, new_g, s, t, NULL);

    Expr* coeff_list = expr_new_function(expr_new_symbol(SYM_List), coeff_args, count);
    free(coeff_args);

    Expr* pair[2] = { g_expr, coeff_list };
    return expr_new_function(expr_new_symbol(SYM_List), pair, 2);
}

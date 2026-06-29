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
    if (!arith_warnings_muted()) {
        char* arg_str = expr_to_string(arg);
        char* call_str = expr_to_string(res);
        fprintf(stderr,
                "ExtendedGCD::exact: Argument %s in %s is not an exact number.\n",
                arg_str ? arg_str : "?", call_str ? call_str : "?");
        free(arg_str);
        free(call_str);
    }
    return NULL;
}

/* ExtendedGCD::egcd -- every argument is exact but at least one is a
 * non-integer rational. */
static Expr* egcd_emit_egcd(Expr* res) {
    if (!arith_warnings_muted()) {
        char* call_str = expr_to_string(res);
        fprintf(stderr,
                "ExtendedGCD::egcd: Arguments in %s should be integers.\n",
                call_str ? call_str : "?");
        free(call_str);
    }
    return NULL;
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
    if (!all_integer) {
        /* All-exact-but-non-integer -> diagnose; anything symbolic stays
         * unevaluated silently so rules / symbolic flow still apply. */
        if (all_exact) return egcd_emit_egcd(res);
        return NULL;
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

/* nt_util.c -- shared GCD/LCM/rational helpers (single_arg_abs_or_copy, rational_like_*, mpz_pair_to_rational_expr, lcm_checked).
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

/* Single-argument GCD/LCM: mathematically the GCD or LCM of one number
 * is its absolute value (GCD[-5] == 5, GCD[-1/2] == 1/2).  Symbolic
 * input is returned unchanged so user-supplied rules / OneIdentity
 * pattern matching still apply (GCD[x] == x). */
Expr* single_arg_abs_or_copy(Expr* arg) {
    int64_t n, d;
    if (is_rational(arg, &n, &d)) {
        return make_rational(llabs(n), llabs(d));
    }
    if (arg->type == EXPR_BIGINT) {
        mpz_t a;
        mpz_init(a);
        mpz_abs(a, arg->data.bigint);
        Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(a));
        mpz_clear(a);
        return r;
    }
    return expr_copy(arg);
}

/* Extract an integer-like or rational-like expression into two
 * initialized mpz_t variables (num, den). Caller owns both; this helper
 * does NOT call mpz_init on its arguments — they must already be init'd.
 * Used by GCD/LCM to support Rational[BigInt, _] without overflowing the
 * int64 rational fold. */
void rational_like_to_mpz_pair(const Expr* e, mpz_t num, mpz_t den) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) {
        expr_to_mpz(e, num);
        mpz_set_ui(den, 1);
        return;
    }
    /* is_rational_like guaranteed: Rational[Integer|BigInt, Integer|BigInt]. */
    expr_to_mpz(e->data.function.args[0], num);
    expr_to_mpz(e->data.function.args[1], den);
}

/* Returns true when at least one component of a rational-like arg is a
 * BigInt (the bare BigInt case, or either part of a Rational). The int64
 * rational fold below cannot handle these; the mpz fold must. */
bool rational_like_needs_bigint(const Expr* e) {
    if (e->type == EXPR_BIGINT) return true;
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Rational &&
        e->data.function.arg_count == 2) {
        return e->data.function.args[0]->type == EXPR_BIGINT ||
               e->data.function.args[1]->type == EXPR_BIGINT;
    }
    return false;
}

/* Build a normalised Integer or Rational expression from a num/den pair
 * supplied as mpz_t. Reduces the pair by their GCD first; if the
 * resulting denominator is 1, returns Integer/BigInt; otherwise builds
 * Rational[num, den]. Inputs are read-only and not freed. */
Expr* mpz_pair_to_rational_expr(const mpz_t num_in, const mpz_t den_in) {
    mpz_t n, d, g;
    mpz_init_set(n, num_in);
    mpz_init_set(d, den_in);
    mpz_init(g);
    mpz_gcd(g, n, d);
    if (mpz_cmp_ui(g, 1) > 0) {
        mpz_divexact(n, n, g);
        mpz_divexact(d, d, g);
    }
    Expr* out;
    if (mpz_cmp_ui(d, 1) == 0) {
        out = expr_bigint_normalize(expr_new_bigint_from_mpz(n));
    } else {
        Expr* en = expr_bigint_normalize(expr_new_bigint_from_mpz(n));
        Expr* ed = expr_bigint_normalize(expr_new_bigint_from_mpz(d));
        Expr* args[2] = { en, ed };
        out = expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
    }
    mpz_clears(n, d, g, NULL);
    return out;
}

/* int64 LCM with overflow detection. Returns LCM(|a|, |b|); sets *overflow
 * when the result does not fit in int64_t (the int64 fast paths fall back to
 * the GMP folds below in that case). Callers handle the zero operands. */
int64_t lcm_checked(int64_t a, int64_t b, bool *overflow) {
    if (a == 0 || b == 0) return 0;
    a = llabs(a);
    b = llabs(b);
    int64_t g = gcd(a, b);
    __int128 prod = (__int128)(a / g) * (__int128)b;
    if (prod > (__int128)INT64_MAX) {
        *overflow = true;
        return 0;
    }
    return (int64_t)prod;
}

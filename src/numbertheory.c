/* numbertheory.c — number-theoretic builtins split out of arithmetic.c.
 *
 * Owns GCD, LCM, ExtendedGCD, PowerMod (and its modular-root helpers),
 * Factorial, Factorial2, FactorialPower, Binomial, PrimitiveRoot,
 * PrimitiveRootList, MultiplicativeOrder and their private helpers.  The
 * core rational/complex constructors, the shared int64 gcd()/lcm() helpers,
 * and the shared numeric predicates (is_infinity_sym, expr_numeric_sign,
 * ...) remain in arithmetic.c (declared in arithmetic.h, which this file
 * includes). */

#include "numbertheory.h"
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
static Expr* single_arg_abs_or_copy(Expr* arg) {
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
static void rational_like_to_mpz_pair(const Expr* e, mpz_t num, mpz_t den) {
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
static bool rational_like_needs_bigint(const Expr* e) {
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
static Expr* mpz_pair_to_rational_expr(const mpz_t num_in, const mpz_t den_in) {
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

    int64_t running_n = 0;
    int64_t running_d = 1;

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
            running_d = lcm(running_d, d);
        }
    }

    return make_rational(running_n, running_d);
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

    int64_t running_n = 0;
    int64_t running_d = 0;

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
                running_n = lcm(running_n, n);
                running_d = gcd(running_d, d);
            }
        }
    }

    return make_rational(running_n, running_d);
}

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





/* ===== PowerMod modular-root helpers ============================
 *
 * The integer-exponent case is handled directly by GMP's mpz_powm /
 * mpz_invert.  The rational-exponent case PowerMod[a, p/q, m] asks for
 * x such that x^q  ≡  a^p  (mod m); GMP has no primitive for modular
 * r-th roots, so we build one:
 *
 *   1. Brute-force scan for tiny moduli (m ≤ MODROOT_BRUTE_LIMIT).
 *   2. Otherwise factor m via internal_factorinteger.
 *   3. For each prime power p^e:
 *        a. solve x_0^r ≡ c (mod p) using
 *             - Tonelli-Shanks when r = 2,
 *             - the closed form  x = c^(r^-1 mod (p-1))  when
 *               gcd(r, p-1) = 1,
 *             - brute force for small primes,
 *        b. Hensel-lift x_0 to mod p^e.
 *   4. Combine the per-prime-power roots with CRT.
 *
 * Returns 0 ("no solution found / unsupported case") rather than failing
 * the whole PowerMod, so the surface expression stays unevaluated and the
 * user sees PowerMod[...] echoed back.
 */

#define MODROOT_BRUTE_LIMIT 1000000UL

static int modroot_brute(mpz_t x, const mpz_t c, const mpz_t r, const mpz_t m) {
    if (!mpz_fits_ulong_p(m)) return 0;
    unsigned long m_ul = mpz_get_ui(m);
    if (m_ul == 0 || m_ul > MODROOT_BRUTE_LIMIT) return 0;
    mpz_t c_mod, ii, ip;
    mpz_inits(c_mod, ii, ip, NULL);
    mpz_mod(c_mod, c, m);
    int found = 0;
    for (unsigned long k = 0; k < m_ul; k++) {
        mpz_set_ui(ii, k);
        mpz_powm(ip, ii, r, m);
        if (mpz_cmp(ip, c_mod) == 0) {
            mpz_set(x, ii);
            found = 1;
            break;
        }
    }
    mpz_clears(c_mod, ii, ip, NULL);
    return found;
}

/* Tonelli-Shanks: x with x^2 ≡ a (mod p), p odd prime.  a assumed already
 * reduced into [0, p).  Returns 0 if a is a quadratic non-residue. */
static int tonelli_shanks(mpz_t x, const mpz_t a, const mpz_t p) {
    if (mpz_sgn(a) == 0) { mpz_set_ui(x, 0); return 1; }

    mpz_t pm1, half, leg;
    mpz_inits(pm1, half, leg, NULL);
    mpz_sub_ui(pm1, p, 1);
    mpz_fdiv_q_2exp(half, pm1, 1);
    mpz_powm(leg, a, half, p);
    if (mpz_cmp_ui(leg, 1) != 0) {
        mpz_clears(pm1, half, leg, NULL);
        return 0;
    }

    /* p ≡ 3 (mod 4) shortcut: x = a^((p+1)/4) mod p. */
    if (mpz_tstbit(p, 1) == 1) {
        mpz_t exp;
        mpz_init(exp);
        mpz_add_ui(exp, p, 1);
        mpz_fdiv_q_2exp(exp, exp, 2);
        mpz_powm(x, a, exp, p);
        mpz_clears(exp, pm1, half, leg, NULL);
        return 1;
    }

    /* Write p - 1 = q * 2^s with q odd. */
    mpz_t q;
    mpz_init_set(q, pm1);
    unsigned long s = 0;
    while (mpz_even_p(q)) { mpz_fdiv_q_2exp(q, q, 1); s++; }

    /* Find a quadratic non-residue z by trial. */
    mpz_t z, leg_z;
    mpz_inits(z, leg_z, NULL);
    mpz_set_ui(z, 2);
    for (;;) {
        mpz_powm(leg_z, z, half, p);
        if (mpz_cmp(leg_z, pm1) == 0) break;
        mpz_add_ui(z, z, 1);
        if (mpz_cmp(z, p) >= 0) {
            mpz_clears(z, leg_z, q, pm1, half, leg, NULL);
            return 0;
        }
    }

    mpz_t M, c, t, R, b, tt, e0;
    mpz_inits(M, c, t, R, b, tt, e0, NULL);
    mpz_set_ui(M, s);
    mpz_powm(c, z, q, p);
    mpz_powm(t, a, q, p);
    mpz_add_ui(e0, q, 1);
    mpz_fdiv_q_2exp(e0, e0, 1);
    mpz_powm(R, a, e0, p);

    int rc = 1;
    while (mpz_cmp_ui(t, 1) != 0) {
        unsigned long i = 0;
        mpz_set(tt, t);
        while (mpz_cmp_ui(tt, 1) != 0) {
            mpz_mul(tt, tt, tt);
            mpz_mod(tt, tt, p);
            i++;
            if (mpz_cmp_ui(M, i) <= 0) { rc = 0; goto ts_done; }
        }
        unsigned long m_ul = mpz_get_ui(M);
        mpz_set(b, c);
        for (unsigned long k = 0; k < m_ul - i - 1; k++) {
            mpz_mul(b, b, b);
            mpz_mod(b, b, p);
        }
        mpz_set_ui(M, i);
        mpz_mul(c, b, b); mpz_mod(c, c, p);
        mpz_mul(t, t, c); mpz_mod(t, t, p);
        mpz_mul(R, R, b); mpz_mod(R, R, p);
    }
    mpz_set(x, R);
ts_done:
    mpz_clears(M, c, t, R, b, tt, e0, z, leg_z, q, pm1, half, leg, NULL);
    return rc;
}

/* x with x^r ≡ a (mod p), p prime.  Returns 0 if no root or unsupported. */
static int rth_root_mod_prime(mpz_t x, const mpz_t a, const mpz_t r,
                              const mpz_t p) {
    mpz_t a_mod;
    mpz_init(a_mod);
    mpz_mod(a_mod, a, p);
    if (mpz_sgn(a_mod) == 0) { mpz_set_ui(x, 0); mpz_clear(a_mod); return 1; }
    if (mpz_cmp_ui(r, 1) == 0) { mpz_set(x, a_mod); mpz_clear(a_mod); return 1; }

    if (mpz_cmp_ui(p, 2) == 0) {
        /* a_mod is 1 here (0 case handled above); 1^r = 1 mod 2. */
        mpz_set_ui(x, 1);
        mpz_clear(a_mod);
        return 1;
    }

    if (mpz_cmp_ui(r, 2) == 0) {
        int rc = tonelli_shanks(x, a_mod, p);
        mpz_clear(a_mod);
        return rc;
    }

    /* General r > 2: try coprime closed form. */
    mpz_t phi, g, r_inv;
    mpz_inits(phi, g, r_inv, NULL);
    mpz_sub_ui(phi, p, 1);
    mpz_gcd(g, r, phi);
    if (mpz_cmp_ui(g, 1) == 0) {
        mpz_invert(r_inv, r, phi);
        mpz_powm(x, a_mod, r_inv, p);
        mpz_clears(phi, g, r_inv, a_mod, NULL);
        return 1;
    }

    /* gcd > 1: need a^(phi/g) ≡ 1 (mod p) for solvability.  When that
     * holds the full Adleman-Manders-Miller algorithm is needed; we fall
     * back to brute force for small primes and otherwise give up. */
    mpz_t phi_g, leg;
    mpz_inits(phi_g, leg, NULL);
    mpz_fdiv_q(phi_g, phi, g);
    mpz_powm(leg, a_mod, phi_g, p);
    int solvable = (mpz_cmp_ui(leg, 1) == 0);
    mpz_clears(phi_g, leg, NULL);

    int rc = 0;
    if (solvable) {
        rc = modroot_brute(x, a_mod, r, p);
    }
    mpz_clears(phi, g, r_inv, a_mod, NULL);
    return rc;
}

/* Newton/Hensel lift: given x0 with x0^r ≡ c (mod p) and gcd(p, r*x0) = 1,
 * compute x with x^r ≡ c (mod p^e). */
static int hensel_lift(mpz_t x, const mpz_t x0, const mpz_t r, const mpz_t c,
                       const mpz_t p, unsigned long e) {
    if (e == 1) { mpz_set(x, x0); return 1; }
    mpz_t xk, fx, fpx, inv, prod, rm1, pk_next;
    mpz_inits(xk, fx, fpx, inv, prod, rm1, pk_next, NULL);
    mpz_set(xk, x0);
    mpz_sub_ui(rm1, r, 1);
    unsigned long k = 1;
    int rc = 1;
    while (k < e) {
        unsigned long next_k = 2 * k;
        if (next_k > e) next_k = e;
        mpz_pow_ui(pk_next, p, next_k);
        mpz_powm(fx, xk, r, pk_next);
        mpz_sub(fx, fx, c);
        mpz_mod(fx, fx, pk_next);
        mpz_powm(fpx, xk, rm1, pk_next);
        mpz_mul(fpx, fpx, r);
        mpz_mod(fpx, fpx, pk_next);
        if (mpz_invert(inv, fpx, pk_next) == 0) { rc = 0; break; }
        mpz_mul(prod, fx, inv);
        mpz_sub(xk, xk, prod);
        mpz_mod(xk, xk, pk_next);
        k = next_k;
    }
    if (rc) mpz_set(x, xk);
    mpz_clears(xk, fx, fpx, inv, prod, rm1, pk_next, NULL);
    return rc;
}

/* x with x^r ≡ c (mod p^e); also write p^e into pe_out. */
static int rth_root_mod_pe(mpz_t x, mpz_t pe_out, const mpz_t c,
                           const mpz_t r, const mpz_t p, unsigned long e) {
    mpz_pow_ui(pe_out, p, e);

    if (modroot_brute(x, c, r, pe_out)) return 1;

    /* Solve mod p first, then Hensel-lift to mod p^e. */
    mpz_t c_p, x0;
    mpz_inits(c_p, x0, NULL);
    mpz_mod(c_p, c, p);
    int rc = rth_root_mod_prime(x0, c_p, r, p);
    if (!rc) { mpz_clears(c_p, x0, NULL); return 0; }

    /* If x0 ≡ 0 mod p but c ≢ 0 mod p^e the lift is the p|c branch we
     * don't implement here.  brute force already covered small pe. */
    if (mpz_sgn(x0) == 0) {
        mpz_t c_mod;
        mpz_init(c_mod);
        mpz_mod(c_mod, c, pe_out);
        if (mpz_sgn(c_mod) == 0) { mpz_set_ui(x, 0); rc = 1; }
        else rc = 0;
        mpz_clears(c_mod, c_p, x0, NULL);
        return rc;
    }
    /* Hensel needs p ∤ r * x0^(r-1).  x0 != 0 mod p; require p ∤ r too. */
    if (mpz_divisible_p(r, p)) { mpz_clears(c_p, x0, NULL); return 0; }

    rc = hensel_lift(x, x0, r, c, p, e);
    mpz_clears(c_p, x0, NULL);
    return rc;
}

/* Top-level: find x in [0, m) with x^r ≡ c (mod m). */
static int modular_root(mpz_t x, const mpz_t c, const mpz_t r, const mpz_t m) {
    if (mpz_cmp_ui(m, 1) == 0) { mpz_set_ui(x, 0); return 1; }
    if (mpz_cmp_ui(m, 1) < 0)  return 0;
    if (mpz_cmp_ui(r, 0) <= 0) return 0;
    if (mpz_cmp_ui(r, 1) == 0) { mpz_mod(x, c, m); return 1; }

    if (modroot_brute(x, c, r, m)) return 1;
    if (mpz_probab_prime_p(m, 25)) return rth_root_mod_prime(x, c, r, m);

    /* Composite m: factor via FactorInteger, solve per prime power, CRT.
     * internal_factorinteger takes ownership of the args array contents,
     * so we must NOT free m_expr afterwards. */
    Expr* args[1] = { expr_bigint_normalize(expr_new_bigint_from_mpz(m)) };
    Expr* fact = internal_factorinteger(args, 1);
    if (!fact || fact->type != EXPR_FUNCTION ||
        fact->data.function.head->type != EXPR_SYMBOL ||
        fact->data.function.head->data.symbol != SYM_List) {
        if (fact) expr_free(fact);
        return 0;
    }

    size_t n = fact->data.function.arg_count;
    mpz_t x_acc, m_acc;
    mpz_inits(x_acc, m_acc, NULL);
    mpz_set_ui(x_acc, 0);
    mpz_set_ui(m_acc, 1);

    int rc = 1;
    for (size_t i = 0; i < n; i++) {
        Expr* pair = fact->data.function.args[i];
        if (pair->type != EXPR_FUNCTION || pair->data.function.arg_count != 2 ||
            !expr_is_integer_like(pair->data.function.args[0]) ||
            pair->data.function.args[1]->type != EXPR_INTEGER ||
            pair->data.function.args[1]->data.integer <= 0) {
            rc = 0; break;
        }
        mpz_t pi, pe, xi;
        mpz_inits(pe, xi, NULL);
        expr_to_mpz(pair->data.function.args[0], pi);
        unsigned long ei = (unsigned long)pair->data.function.args[1]->data.integer;
        if (!rth_root_mod_pe(xi, pe, c, r, pi, ei)) {
            mpz_clears(pi, pe, xi, NULL);
            rc = 0; break;
        }
        if (i == 0) {
            mpz_set(x_acc, xi);
            mpz_set(m_acc, pe);
        } else {
            mpz_t inv, diff, term;
            mpz_inits(inv, diff, term, NULL);
            if (mpz_invert(inv, m_acc, pe) == 0) {
                mpz_clears(inv, diff, term, pi, pe, xi, NULL);
                rc = 0; break;
            }
            mpz_sub(diff, xi, x_acc);
            mpz_mul(term, diff, inv);
            mpz_mod(term, term, pe);
            mpz_mul(term, term, m_acc);
            mpz_add(x_acc, x_acc, term);
            mpz_mul(m_acc, m_acc, pe);
            mpz_mod(x_acc, x_acc, m_acc);
            mpz_clears(inv, diff, term, NULL);
        }
        mpz_clears(pi, pe, xi, NULL);
    }
    expr_free(fact);
    if (rc) mpz_set(x, x_acc);
    mpz_clears(x_acc, m_acc, NULL);
    return rc;
}

/* PowerMod[a, b, m].  b may be an integer (incl. negative for the modular
 * inverse path) or Rational[p,q] (modular q-th root of a^p mod m). */
Expr* builtin_powermod(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;

    Expr* a_expr = res->data.function.args[0];
    Expr* b_expr = res->data.function.args[1];
    Expr* m_expr = res->data.function.args[2];

    if (!expr_is_integer_like(a_expr) || !expr_is_integer_like(m_expr)) return NULL;

    mpz_t a, m;
    expr_to_mpz(a_expr, a);
    expr_to_mpz(m_expr, m);
    if (mpz_sgn(m) == 0) { mpz_clears(a, m, NULL); return NULL; }
    mpz_abs(m, m);

    /* --- Integer (or BigInt) exponent: GMP fast path. ----------------- */
    if (expr_is_integer_like(b_expr)) {
        mpz_t b, r;
        expr_to_mpz(b_expr, b);
        mpz_init(r);
        if (mpz_sgn(b) < 0) {
            mpz_t inv, b_neg;
            mpz_inits(inv, b_neg, NULL);
            if (mpz_invert(inv, a, m) == 0) {
                mpz_clears(inv, b_neg, b, a, m, r, NULL);
                return NULL;
            }
            mpz_neg(b_neg, b);
            mpz_powm(r, inv, b_neg, m);
            mpz_clears(inv, b_neg, NULL);
        } else {
            mpz_powm(r, a, b, m);
        }
        Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
        mpz_clears(b, a, m, r, NULL);
        return out;
    }

    /* --- Rational exponent p/q: modular q-th root of a^p. ------------- */
    if (b_expr->type == EXPR_FUNCTION && b_expr->data.function.head &&
        b_expr->data.function.head->type == EXPR_SYMBOL &&
        b_expr->data.function.head->data.symbol == SYM_Rational &&
        b_expr->data.function.arg_count == 2 &&
        expr_is_integer_like(b_expr->data.function.args[0]) &&
        expr_is_integer_like(b_expr->data.function.args[1])) {
        mpz_t p, q, c, root;
        expr_to_mpz(b_expr->data.function.args[0], p);
        expr_to_mpz(b_expr->data.function.args[1], q);
        mpz_inits(c, root, NULL);

        if (mpz_sgn(q) <= 0) {
            mpz_clears(p, q, c, root, a, m, NULL);
            return NULL;
        }

        /* c = a^p mod m, handling negative p via modular inverse of a. */
        if (mpz_sgn(p) < 0) {
            mpz_t inv, p_neg;
            mpz_inits(inv, p_neg, NULL);
            if (mpz_invert(inv, a, m) == 0) {
                mpz_clears(inv, p_neg, p, q, c, root, a, m, NULL);
                return NULL;
            }
            mpz_neg(p_neg, p);
            mpz_powm(c, inv, p_neg, m);
            mpz_clears(inv, p_neg, NULL);
        } else {
            mpz_powm(c, a, p, m);
        }

        if (!modular_root(root, c, q, m)) {
            mpz_clears(p, q, c, root, a, m, NULL);
            return NULL;
        }
        Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(root));
        mpz_clears(p, q, c, root, a, m, NULL);
        return out;
    }

    mpz_clears(a, m, NULL);
    return NULL;
}

Expr* builtin_factorial(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* Machine Real: Factorial[x] = Gamma[x + 1] via libm tgamma. */
    if (arg->type == EXPR_REAL) {
        double v = arg->data.real;
        return expr_new_real(tgamma(v + 1.0));
    }
#ifdef USE_MPFR
    /* MPFR Real: same identity at full input precision. */
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        mpfr_t shifted;
        mpfr_init2(shifted, prec);
        mpfr_add_ui(shifted, arg->data.mpfr, 1, MPFR_RNDN);
        Expr* result = expr_new_mpfr_bits(prec);
        mpfr_gamma(result->data.mpfr, shifted, MPFR_RNDN);
        mpfr_clear(shifted);
        return result;
    }
#endif
    /* BigInt: factorial of a value that exceeds int64 is astronomical
     * (1e20! has ~10^21 digits) and would exhaust memory. Mathematica
     * leaves it symbolic for the same reason. */

    int64_t n, d;
    if (is_rational(arg, &n, &d)) {
        if (d == 1) {
            if (n < 0) return expr_new_symbol(SYM_ComplexInfinity);
            if (n <= 20) {
                int64_t f = 1;
                for (int64_t i = 2; i <= n; i++) f *= i;
                return expr_new_integer(f);
            } else {
                mpz_t result;
                mpz_init(result);
                mpz_fac_ui(result, (unsigned long)n);
                Expr* r = expr_new_bigint_from_mpz(result);
                mpz_clear(result);
                return r;
            }
        } else if (d == 2 || d == -2) {
            if (d == -2) { n = -n; }
            /* Half-integer Factorial: (n/2)! = coeff * Sqrt[Pi], where the
             * coefficient is built in GMP -- the odd double factorial and the
             * 2^k denominator both blow past int64 well before the argument is
             * large (e.g. Gamma[201/2] needs 199!! / 2^100 ~ 10^157). */
            mpz_t num, den;
            mpz_init_set_ui(num, 1);
            mpz_init_set_ui(den, 1);

            if (n > 0) {
                for (int64_t i = n; i >= 1; i -= 2) mpz_mul_si(num, num, (long)i);
                mpz_mul_2exp(den, den, (unsigned long)((n + 1) / 2)); /* 2^((n+1)/2) */
            } else {
                for (int64_t i = n + 2; i <= -1; i += 2) {
                    mpz_mul_ui(num, num, 2);
                    mpz_mul_si(den, den, (long)i);
                }
            }
            /* Keep the denominator positive (negative half-integers can flip it). */
            if (mpz_sgn(den) < 0) { mpz_neg(num, num); mpz_neg(den, den); }

            Expr* coeff = mpz_pair_to_rational_expr(num, den);
            mpz_clears(num, den, NULL);
            if (!coeff) coeff = expr_new_integer(0);

            Expr* pi_sym = expr_new_symbol(SYM_Pi);
            Expr* half = make_rational(1, 2);
            Expr* sqrt_pi = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){pi_sym, half}, 2));
            
            if (coeff->type == EXPR_INTEGER && coeff->data.integer == 1) {
                expr_free(coeff);
                return sqrt_pi;
            } else {
                return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){coeff, sqrt_pi}, 2));
            }
        }
    }

    return NULL;
}

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

/* Factorial2[n] -- the double factorial.
 *
 * For non-negative integer n:
 *   n!! = n * (n-2) * (n-4) * ...    (terminating at 1 if n odd, 2 if n even).
 *   0!! = 1, (-1)!! = 1.
 * For negative odd integers n <= -3:
 *   n!! = ComplexInfinity (poles, like Factorial of negative integers).
 *
 * Symbolic / non-integer arguments return NULL (the call stays as-is)
 * so the surrounding evaluator can keep them in expression form.
 *
 * Implemented in two paths: int64 fast-path for small n (<= 30 -- fits in
 * int64 for all valid inputs since 31!! ~= 7.9e9 fits), and a GMP path
 * for larger n. */
Expr* builtin_factorial2(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* Special-case the small literal -1 / 0 inputs. */
    if (arg->type == EXPR_INTEGER) {
        int64_t n = arg->data.integer;
        if (n == -1 || n == 0) return expr_new_integer(1);
        if (n < 0) {
            /* Even negative integers: ComplexInfinity (pole). Odd negative
             * integers below -1 also diverge. */
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        if (n <= 30) {
            int64_t f = 1;
            for (int64_t i = n; i >= 1; i -= 2) f *= i;
            return expr_new_integer(f);
        }
        /* GMP path. */
        mpz_t r, k;
        mpz_init_set_ui(r, 1);
        mpz_init(k);
        for (int64_t i = n; i >= 1; i -= 2) {
            mpz_set_si(k, (long)i);
            mpz_mul(r, r, k);
        }
        Expr* out = expr_new_bigint_from_mpz(r);
        mpz_clears(r, k, NULL);
        return out;
    }

    if (arg->type == EXPR_BIGINT) {
        /* For arguments large enough to be BigInt, use GMP. We only
         * handle non-negative inputs; negative BigInt is ComplexInfinity. */
        if (mpz_sgn(arg->data.bigint) < 0) {
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        if (mpz_sgn(arg->data.bigint) == 0) return expr_new_integer(1);
        /* Cap at a sane size: if mpz_fits_ulong_p fails the input is too
         * large for our explicit loop. Return NULL to leave it symbolic. */
        if (!mpz_fits_ulong_p(arg->data.bigint)) return NULL;
        unsigned long n = mpz_get_ui(arg->data.bigint);
        mpz_t r, k;
        mpz_init_set_ui(r, 1);
        mpz_init(k);
        for (unsigned long i = n; i >= 1; i -= 2) {
            mpz_set_ui(k, i);
            mpz_mul(r, r, k);
            if (i < 2) break; /* guard underflow when i==1 */
        }
        Expr* out = expr_new_bigint_from_mpz(r);
        mpz_clears(r, k, NULL);
        return out;
    }

    return NULL; /* symbolic -- leave Factorial2[expr] alone */
}

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



/* =====================================================================
 * PrimitiveRoot / PrimitiveRootList
 *
 *   PrimitiveRoot[n]        a primitive root of n
 *   PrimitiveRoot[n, k]     smallest primitive root of n >= k
 *   PrimitiveRootList[n]    sorted list of primitive roots of n in [1, n-1]
 *
 * A primitive root of n is a generator of the multiplicative group
 * (Z/nZ)^*, which is cyclic iff n in {1, 2, 4} or n = p^k or n = 2 p^k
 * for an odd prime p and k >= 1.  For non-cyclic n the call is left
 * unevaluated (PrimitiveRoot) or returns {} (PrimitiveRootList).
 *
 * All arithmetic is on GMP mpz_t so machine and big integers share the
 * same code path; int64 inputs are coerced via expr_to_mpz.
 * ===================================================================*/

/* Naive prime test for unsigned long; only used to skip composite trial
 * exponents in pr_decompose_odd_prime_power.  q <= log2(n) is tiny in
 * practice (a few hundred at most for any realistic input). */
static bool pr_ulong_is_prime(unsigned long q) {
    if (q < 2) return false;
    if (q < 4) return true;
    if ((q & 1UL) == 0) return false;
    for (unsigned long d = 3; d * d <= q; d += 2) {
        if (q % d == 0) return false;
    }
    return true;
}

/* Decompose an odd integer n >= 3 as p^k with p prime.  Returns true on
 * success; on success p_out is set to p and *k_out to k.  Returns false
 * if n is composite but not a prime power, or if the base turns out to
 * be 2 (n was supposed to be odd).
 *
 * Algorithm: iteratively strip prime exponents q from the current value
 * by computing the exact q-th root.  After O(omega(k)) reductions the
 * current value is the prime base p.  Each mpz_root is O(M(size)) so
 * the whole routine is dominated by the primality test on the final p. */
static bool pr_decompose_odd_prime_power(const mpz_t n, mpz_t p_out, uint64_t* k_out) {
    if (mpz_cmp_ui(n, 1) <= 0) return false;
    mpz_t current;
    mpz_init_set(current, n);
    uint64_t total_k = 1;

    while (mpz_probab_prime_p(current, 25) == 0) {
        unsigned long log2_cur = (unsigned long)mpz_sizeinbase(current, 2);
        bool factored = false;
        for (unsigned long q = 2; q <= log2_cur; q++) {
            if (!pr_ulong_is_prime(q)) continue;
            mpz_t cand;
            mpz_init(cand);
            int exact = mpz_root(cand, current, q);
            if (exact) {
                mpz_set(current, cand);
                total_k *= q;
                mpz_clear(cand);
                factored = true;
                break;
            }
            mpz_clear(cand);
        }
        if (!factored) {
            mpz_clear(current);
            return false;
        }
    }
    /* current is prime.  Reject p == 2 since we only recognise odd prime
     * power moduli here; the 2 / 4 cases are handled at the caller. */
    if (mpz_cmp_ui(current, 2) == 0) {
        mpz_clear(current);
        return false;
    }
    mpz_set(p_out, current);
    *k_out = total_k;
    mpz_clear(current);
    return true;
}

/* Classify n as a cyclic modulus.  On success sets:
 *   *is_two_pk_out  true iff n = 2 p^k (k >= 1, odd prime p)
 *   p_out, *k_out   the odd prime base and exponent (only meaningful for
 *                   the n = p^k and n = 2 p^k cases; for n in {2, 4}
 *                   they are not written).
 * Returns one of CYCLIC_NONE / CYCLIC_TWO / CYCLIC_FOUR / CYCLIC_PRIME_POWER
 * / CYCLIC_TWO_PRIME_POWER.  n must be a positive mpz_t. */
typedef enum {
    PR_CYC_NONE = 0,
    PR_CYC_TWO,            /* n == 2 */
    PR_CYC_FOUR,           /* n == 4 */
    PR_CYC_PRIME_POWER,    /* n == p^k, p odd prime, k >= 1 */
    PR_CYC_TWO_PRIME_POWER /* n == 2 p^k */
} pr_cyclic_kind;

static pr_cyclic_kind pr_classify(const mpz_t n, mpz_t p_out, uint64_t* k_out) {
    if (mpz_cmp_ui(n, 2) == 0) return PR_CYC_TWO;
    if (mpz_cmp_ui(n, 4) == 0) return PR_CYC_FOUR;
    if (mpz_cmp_ui(n, 1) <= 0) return PR_CYC_NONE;

    if (mpz_odd_p(n)) {
        if (pr_decompose_odd_prime_power(n, p_out, k_out)) return PR_CYC_PRIME_POWER;
        return PR_CYC_NONE;
    }
    /* even: n must be 2 * (odd prime power).  Strip exactly one factor of 2. */
    if (mpz_divisible_2exp_p(n, 2)) return PR_CYC_NONE; /* 4 | n => non-cyclic */
    mpz_t half;
    mpz_init(half);
    mpz_tdiv_q_2exp(half, n, 1);
    if (mpz_cmp_ui(half, 1) == 0) {
        /* n == 2, already handled above. */
        mpz_clear(half);
        return PR_CYC_NONE;
    }
    bool ok = pr_decompose_odd_prime_power(half, p_out, k_out);
    mpz_clear(half);
    return ok ? PR_CYC_TWO_PRIME_POWER : PR_CYC_NONE;
}

/* Maximum number of distinct prime factors of phi(n) we expect.  phi for
 * any cyclic modulus is bounded by n, and omega(m) <= ~15 for m <= 10^18,
 * <= ~50 for m <= 10^60, <= ~330 for m <= 10^900.  256 leaves headroom
 * for moderately huge p while staying stack-friendly. */
#define PR_MAX_DISTINCT_PRIMES 256

/* Collect distinct prime divisors of m (m >= 1) into out[0..*count).
 * out has capacity PR_MAX_DISTINCT_PRIMES; *count must be 0 on entry and
 * each written out[i] is mpz_init_set by the routine.  Returns true on
 * success.  Caller must mpz_clear out[0..*count) on both true and false.
 *
 * Strategy:
 *   - Trial-divide by primes up to PR_TRIAL_LIMIT (covers p - 1 for any
 *     p up to PR_TRIAL_LIMIT^2 ~ 1e10 in one pass).
 *   - If a non-trivial cofactor remains, descend via Pollard rho until
 *     all factors are prime.
 * This avoids pulling in the full FactorInteger pipeline for what is
 * almost always a small number (p - 1). */
#define PR_TRIAL_LIMIT 100000UL

static void pr_pollard_rho(mpz_t f_out, const mpz_t n_in, unsigned long c) {
    mpz_t x, y, d, diff, n;
    mpz_inits(x, y, d, diff, n, NULL);
    mpz_set(n, n_in);
    mpz_set_ui(x, 2);
    mpz_set_ui(y, 2);
    mpz_set_ui(d, 1);
    unsigned long steps = 0;
    while (mpz_cmp_ui(d, 1) == 0 && steps < 1000000UL) {
        /* x = x^2 + c mod n */
        mpz_mul(x, x, x); mpz_add_ui(x, x, c); mpz_mod(x, x, n);
        /* y = (y^2 + c)^2 + c mod n  (advance two steps) */
        mpz_mul(y, y, y); mpz_add_ui(y, y, c); mpz_mod(y, y, n);
        mpz_mul(y, y, y); mpz_add_ui(y, y, c); mpz_mod(y, y, n);
        mpz_sub(diff, x, y); mpz_abs(diff, diff);
        mpz_gcd(d, diff, n);
        steps++;
    }
    if (mpz_cmp(d, n) == 0 || mpz_cmp_ui(d, 1) == 0) {
        mpz_set_ui(f_out, 0); /* failure signal */
    } else {
        mpz_set(f_out, d);
    }
    mpz_clears(x, y, d, diff, n, NULL);
}

static bool pr_add_unique(mpz_t* out, size_t* count, const mpz_t p) {
    for (size_t i = 0; i < *count; i++) {
        if (mpz_cmp(out[i], p) == 0) return true;
    }
    if (*count >= PR_MAX_DISTINCT_PRIMES) return false;
    mpz_init_set(out[*count], p);
    (*count)++;
    return true;
}

static bool pr_distinct_primes_recursive(mpz_t n, mpz_t* out, size_t* count) {
    if (mpz_cmp_ui(n, 1) <= 0) return true;
    if (mpz_probab_prime_p(n, 25)) return pr_add_unique(out, count, n);
    mpz_t f, q;
    mpz_inits(f, q, NULL);
    bool ok = false;
    for (unsigned long c = 1; c < 64 && !ok; c++) {
        pr_pollard_rho(f, n, c);
        if (mpz_sgn(f) != 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
            ok = true;
        }
    }
    if (!ok) {
        mpz_clears(f, q, NULL);
        return false;
    }
    mpz_divexact(q, n, f);
    bool result = pr_distinct_primes_recursive(f, out, count) &&
                  pr_distinct_primes_recursive(q, out, count);
    mpz_clears(f, q, NULL);
    return result;
}

/* Append distinct prime factors of m_in to out[*count..]; preserves any
 * entries already present (callers chain this with pr_add_unique). */
static bool pr_collect_distinct_primes(const mpz_t m_in, mpz_t* out, size_t* count) {
    if (mpz_cmp_ui(m_in, 1) <= 0) return true;
    mpz_t m;
    mpz_init_set(m, m_in);
    for (unsigned long p = 2; p < PR_TRIAL_LIMIT && mpz_cmp_ui(m, 1) > 0; p++) {
        if (mpz_divisible_ui_p(m, p)) {
            mpz_t pz;
            mpz_init_set_ui(pz, p);
            if (!pr_add_unique(out, count, pz)) {
                mpz_clear(pz);
                mpz_clear(m);
                return false;
            }
            mpz_clear(pz);
            while (mpz_divisible_ui_p(m, p)) mpz_divexact_ui(m, m, p);
        }
    }
    bool ok = true;
    if (mpz_cmp_ui(m, 1) > 0) ok = pr_distinct_primes_recursive(m, out, count);
    mpz_clear(m);
    return ok;
}

/* Build phi(n) and the distinct prime factors of phi(n) for a cyclic
 * modulus described by `kind`, `p`, `k`.  Returns true on success.
 * phi_out is mpz_init'd by the caller; primes[i] are mpz_init_set by
 * this routine and must be cleared by the caller. */
static bool pr_phi_and_primes(pr_cyclic_kind kind, const mpz_t p, uint64_t k,
                              mpz_t phi_out, mpz_t* primes, size_t* nprimes) {
    *nprimes = 0;
    switch (kind) {
        case PR_CYC_TWO:
            mpz_set_ui(phi_out, 1);
            return true;
        case PR_CYC_FOUR: {
            mpz_set_ui(phi_out, 2);
            mpz_t two;
            mpz_init_set_ui(two, 2);
            bool ok = pr_add_unique(primes, nprimes, two);
            mpz_clear(two);
            return ok;
        }
        case PR_CYC_PRIME_POWER:
        case PR_CYC_TWO_PRIME_POWER: {
            /* phi(p^k)   = p^(k-1) (p - 1)
             * phi(2 p^k) = phi(2) phi(p^k) = p^(k-1) (p - 1)  (same value) */
            mpz_t pm1, pkm1;
            mpz_inits(pm1, pkm1, NULL);
            mpz_sub_ui(pm1, p, 1);
            if (k == 1) {
                mpz_set_ui(pkm1, 1);
            } else {
                mpz_pow_ui(pkm1, p, (unsigned long)(k - 1));
            }
            mpz_mul(phi_out, pkm1, pm1);

            /* Distinct prime factors of phi: {p if k >= 2} ∪ primes(p - 1). */
            bool ok = true;
            if (k >= 2) ok = pr_add_unique(primes, nprimes, p);
            if (ok) ok = pr_collect_distinct_primes(pm1, primes, nprimes);
            mpz_clears(pm1, pkm1, NULL);
            return ok;
        }
        default:
            return false;
    }
}

/* True iff g is a primitive root of n, given phi = phi(n) and the
 * distinct prime divisors of phi.  Pre: 1 <= g < n + something
 * (we always reduce modulo n).  Uses mpz_powm for the order test. */
static bool pr_is_primitive_root(const mpz_t g_in, const mpz_t n, const mpz_t phi,
                                 const mpz_t* primes, size_t nprimes) {
    mpz_t g, gcd_gn;
    mpz_inits(g, gcd_gn, NULL);
    mpz_mod(g, g_in, n);
    if (mpz_sgn(g) == 0) {
        mpz_clears(g, gcd_gn, NULL);
        return false;
    }
    mpz_gcd(gcd_gn, g, n);
    if (mpz_cmp_ui(gcd_gn, 1) != 0) {
        mpz_clears(g, gcd_gn, NULL);
        return false;
    }
    /* Phi == 1 (n in {1, 2}): only the residue 1 generates the trivial
     * group, so any unit is a "primitive root".  All units of Z/2Z is just
     * {1}, which we've already accepted by gcd above. */
    if (mpz_cmp_ui(phi, 1) == 0) {
        mpz_clears(g, gcd_gn, NULL);
        return true;
    }
    mpz_t exp, power;
    mpz_inits(exp, power, NULL);
    bool is_pr = true;
    for (size_t i = 0; i < nprimes; i++) {
        mpz_divexact(exp, phi, primes[i]);
        mpz_powm(power, g, exp, n);
        if (mpz_cmp_ui(power, 1) == 0) { is_pr = false; break; }
    }
    mpz_clears(g, gcd_gn, exp, power, NULL);
    return is_pr;
}

/* Find the smallest primitive root of n that is >= start.  On success
 * writes the root to g_out and returns true.  Caller pre-inits g_out.
 *
 * The density of primitive roots in [1, n-1] is phi(phi(n)) / phi(n),
 * which for typical p is on the order of 1 / ln(ln(p)), so the search
 * terminates after O(ln(ln(p))) candidates on average. */
static bool pr_smallest_primitive_root(mpz_t g_out, const mpz_t n, const mpz_t phi,
                                       const mpz_t* primes, size_t nprimes,
                                       const mpz_t start) {
    mpz_t g;
    mpz_init(g);
    if (mpz_sgn(start) < 1) mpz_set_ui(g, 1); else mpz_set(g, start);

    /* For tiny n (2, 4) the loop terminates immediately. */
    while (mpz_cmp(g, n) < 0 || mpz_cmp_ui(n, 1) <= 0) {
        if (pr_is_primitive_root(g, n, phi, primes, nprimes)) {
            mpz_set(g_out, g);
            mpz_clear(g);
            return true;
        }
        mpz_add_ui(g, g, 1);
        /* Defensive cap: a primitive root always exists in [1, n-1] for
         * cyclic n, so we should never iterate beyond n.  The break below
         * catches the n == 2 case (where the only PR is 1 and start >= 2
         * would otherwise loop forever). */
        if (mpz_cmp(g, n) >= 0 && mpz_cmp(start, n) >= 0) break;
    }
    mpz_clear(g);
    return false;
}

/* Emit `PrimitiveRoot::argt: PrimitiveRoot called with N arguments; 1 or
 * 2 arguments are expected.` to stderr and return NULL. */
static Expr* pr_emit_argt(size_t argc) {
    fprintf(stderr,
            "PrimitiveRoot::argt: PrimitiveRoot called with %zu argument%s; "
            "1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Emit `PrimitiveRoot::intg: Integer greater than 1 expected at position
 * <pos> in <call>.` */
static Expr* pr_emit_intg(size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "PrimitiveRoot::intg: Integer greater than 1 expected at position "
            "%zu in %s.\n",
            pos, call_str ? call_str : "?");
    free(call_str);
    return NULL;
}

/* Emit `PrimitiveRootList::argx: PrimitiveRootList called with N arguments;
 * 1 argument is expected.` */
static Expr* prl_emit_argx(size_t argc) {
    fprintf(stderr,
            "PrimitiveRootList::argx: PrimitiveRootList called with %zu argument%s; "
            "1 argument is expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_primitiveroot(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return pr_emit_argt(argc);

    Expr* n_expr = res->data.function.args[0];
    if (!expr_is_integer_like(n_expr)) {
        /* Numeric non-integer (Real, Rational, Complex, ...) triggers the
         * Mathematica-style ::intg diagnostic; symbolic args flow through
         * silently so user-supplied DownValues / pattern matching apply. */
        if (expr_is_numeric_like(n_expr)) return pr_emit_intg(1, res);
        return NULL;
    }

    mpz_t n;
    expr_to_mpz(n_expr, n);
    if (mpz_cmp_ui(n, 2) < 0) {
        mpz_clear(n);
        return pr_emit_intg(1, res);
    }

    mpz_t p; uint64_t k;
    mpz_init(p);
    pr_cyclic_kind kind = pr_classify(n, p, &k);
    if (kind == PR_CYC_NONE) {
        mpz_clears(n, p, NULL);
        return NULL;
    }

    mpz_t phi;
    mpz_init(phi);
    mpz_t primes[PR_MAX_DISTINCT_PRIMES];
    size_t nprimes = 0;
    if (!pr_phi_and_primes(kind, p, k, phi, primes, &nprimes)) {
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        mpz_clears(n, p, phi, NULL);
        return NULL;
    }

    /* For 2-arg PrimitiveRoot[n, k]: forward search from k.
     * For 1-arg PrimitiveRoot[n]:
     *   - n = 2          -> 1
     *   - n = 4          -> 3
     *   - n = p^k        -> smallest primitive root of n
     *   - n = 2 p^k      -> g if odd else g + p^k, where g is the
     *                       smallest primitive root of p^k
     */
    Expr* out = NULL;
    if (argc == 2) {
        Expr* k_expr = res->data.function.args[1];
        if (!expr_is_integer_like(k_expr)) {
            if (expr_is_numeric_like(k_expr)) {
                for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
                mpz_clears(n, p, phi, NULL);
                return pr_emit_intg(2, res);
            }
            for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
            mpz_clears(n, p, phi, NULL);
            return NULL;
        }
        mpz_t start, g;
        mpz_inits(start, g, NULL);
        expr_to_mpz(k_expr, start);
        if (mpz_sgn(start) < 0) mpz_set_ui(start, 0);
        bool found = pr_smallest_primitive_root(g, n, phi, primes, nprimes, start);
        if (found) out = expr_bigint_normalize(expr_new_bigint_from_mpz(g));
        mpz_clears(start, g, NULL);
    } else {
        switch (kind) {
            case PR_CYC_TWO:
                out = expr_new_integer(1);
                break;
            case PR_CYC_FOUR:
                out = expr_new_integer(3);
                break;
            case PR_CYC_PRIME_POWER: {
                mpz_t g, two;
                mpz_inits(g, two, NULL);
                mpz_set_ui(two, 2);
                bool found = pr_smallest_primitive_root(g, n, phi, primes, nprimes, two);
                if (found) out = expr_bigint_normalize(expr_new_bigint_from_mpz(g));
                mpz_clears(g, two, NULL);
                break;
            }
            case PR_CYC_TWO_PRIME_POWER: {
                /* Smallest PR of p^k.  If odd, return as-is; else add p^k. */
                mpz_t pk, g, two;
                mpz_inits(pk, g, two, NULL);
                mpz_pow_ui(pk, p, (unsigned long)k);
                mpz_set_ui(two, 2);
                bool found = pr_smallest_primitive_root(g, pk, phi, primes, nprimes, two);
                if (found) {
                    if (mpz_even_p(g)) mpz_add(g, g, pk);
                    out = expr_bigint_normalize(expr_new_bigint_from_mpz(g));
                }
                mpz_clears(pk, g, two, NULL);
                break;
            }
            default:
                break;
        }
    }

    for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
    mpz_clears(n, p, phi, NULL);
    return out;
}

/* Comparator for sorting an mpz_t* array via qsort. */
static int pr_mpz_cmp_qsort(const void* a, const void* b) {
    const mpz_t* ma = (const mpz_t*)a;
    const mpz_t* mb = (const mpz_t*)b;
    return mpz_cmp(*ma, *mb);
}

Expr* builtin_primitiverootlist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return prl_emit_argx(argc);

    Expr* n_expr = res->data.function.args[0];
    if (!expr_is_integer_like(n_expr)) {
        /* Per spec: PrimitiveRootList[non-integer] stays unevaluated with
         * no diagnostic (Real / Complex / Symbolic all flow through). */
        return NULL;
    }

    mpz_t n;
    expr_to_mpz(n_expr, n);
    if (mpz_cmp_ui(n, 1) <= 0) {
        mpz_clear(n);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    mpz_t p; uint64_t k;
    mpz_init(p);
    pr_cyclic_kind kind = pr_classify(n, p, &k);
    if (kind == PR_CYC_NONE) {
        mpz_clears(n, p, NULL);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    mpz_t phi;
    mpz_init(phi);
    mpz_t primes[PR_MAX_DISTINCT_PRIMES];
    size_t nprimes = 0;
    if (!pr_phi_and_primes(kind, p, k, phi, primes, &nprimes)) {
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        mpz_clears(n, p, phi, NULL);
        return NULL;
    }

    /* The enumeration is g^i mod n for i in [1, phi] coprime to phi;
     * the number of such i is phi(phi(n)).  Guard against pathological
     * sizes: if phi doesn't fit in an unsigned long we cannot index the
     * loop, and even when it does the resulting list could be huge. */
    if (!mpz_fits_ulong_p(phi)) {
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        mpz_clears(n, p, phi, NULL);
        return NULL;
    }
    unsigned long phi_ul = mpz_get_ui(phi);

    /* Smallest primitive root of n itself (not of p^k) — for the n = 2 p^k
     * case we want the residue in [1, n-1], which lets us walk g^i mod n
     * directly. */
    mpz_t g, two;
    mpz_inits(g, two, NULL);
    mpz_set_ui(two, 1);
    bool have_g = pr_smallest_primitive_root(g, n, phi, primes, nprimes, two);
    if (!have_g) {
        mpz_clears(g, two, NULL);
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        mpz_clears(n, p, phi, NULL);
        return NULL;
    }

    /* Walk g^i mod n iteratively; collect i with gcd(i, phi) == 1. */
    mpz_t cur, i_mpz, g_local;
    mpz_inits(cur, i_mpz, g_local, NULL);
    mpz_set_ui(cur, 1);
    mpz_set(g_local, g);

    mpz_t* roots = (mpz_t*)malloc(sizeof(mpz_t) * (phi_ul + 1));
    if (!roots) {
        mpz_clears(g, two, cur, i_mpz, g_local, NULL);
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        mpz_clears(n, p, phi, NULL);
        return NULL;
    }
    size_t root_count = 0;

    for (unsigned long i = 1; i <= phi_ul; i++) {
        mpz_mul(cur, cur, g_local);
        mpz_mod(cur, cur, n);
        unsigned long ggcd;
        if (phi_ul == 0) {
            ggcd = i;
        } else {
            /* gcd(i, phi_ul) in unsigned long. */
            unsigned long a = i, b = phi_ul;
            while (b) { unsigned long t = a % b; a = b; b = t; }
            ggcd = a;
        }
        if (ggcd == 1) {
            mpz_init_set(roots[root_count], cur);
            root_count++;
        }
    }

    qsort(roots, root_count, sizeof(mpz_t), pr_mpz_cmp_qsort);

    Expr** items = NULL;
    if (root_count > 0) {
        items = (Expr**)malloc(sizeof(Expr*) * root_count);
        for (size_t i = 0; i < root_count; i++) {
            items[i] = expr_bigint_normalize(expr_new_bigint_from_mpz(roots[i]));
        }
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, root_count);
    free(items);

    for (size_t i = 0; i < root_count; i++) mpz_clear(roots[i]);
    free(roots);

    mpz_clears(g, two, cur, i_mpz, g_local, NULL);
    for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
    mpz_clears(n, p, phi, NULL);
    return out;
}

/* =====================================================================
 * MultiplicativeOrder
 *
 *   MultiplicativeOrder[k, n]
 *       smallest positive m with k^m == 1 (mod n).
 *   MultiplicativeOrder[k, n, {r1, r2, ...}]
 *       smallest positive m with k^m == r_i (mod n) for some i.
 *
 * Returns unevaluated when:
 *   - n is 0 (no group);
 *   - gcd(k, n) != 1 (k is not a unit mod n, no finite order exists);
 *   - in the 3-arg form, no power of k lands in the residue set, or the
 *     order is too large to iterate over.
 * All arithmetic is on GMP mpz_t so bignum k and n share the same code
 * path with machine integers.  Negative n is treated as |n|; negative or
 * out-of-range k is reduced modulo n.
 * ===================================================================*/

/* phi(n) for n >= 1.  Uses the pr_collect_distinct_primes helper above to
 * gather the distinct prime factors of n, then applies the multiplicative
 * identity phi(n) = n * prod_{p | n} (1 - 1/p), iteratively realised as
 * (phi / p) * (p - 1).  This avoids needing the exact multiplicities. */
static bool mo_eulerphi_mpz(const mpz_t n, mpz_t phi_out) {
    if (mpz_sgn(n) <= 0) return false;
    if (mpz_cmp_ui(n, 1) == 0) { mpz_set_ui(phi_out, 1); return true; }
    mpz_t primes[PR_MAX_DISTINCT_PRIMES];
    size_t nprimes = 0;
    if (!pr_collect_distinct_primes(n, primes, &nprimes)) {
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        return false;
    }
    mpz_set(phi_out, n);
    mpz_t pm1;
    mpz_init(pm1);
    for (size_t i = 0; i < nprimes; i++) {
        mpz_divexact(phi_out, phi_out, primes[i]);
        mpz_sub_ui(pm1, primes[i], 1);
        mpz_mul(phi_out, phi_out, pm1);
    }
    mpz_clear(pm1);
    for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
    return true;
}

/* Multiplicative order of k_in modulo n_in.  Pre: order_out is mpz_init'd
 * by the caller; n_in / k_in are read-only.  Returns true on success.
 * Returns false (without setting order_out) when n is zero, k is not a
 * unit modulo n, or phi(n) cannot be factored.
 *
 * Strategy: any element order divides phi(n).  Compute phi, then for each
 * distinct prime p | phi successively divide phi by p as long as the
 * residual exponent still maps k to 1.  After the loop, the residual is
 * exactly the order. */
static bool mo_order_mpz(const mpz_t k_in, const mpz_t n_in, mpz_t order_out) {
    if (mpz_sgn(n_in) == 0) return false;
    mpz_t n;
    mpz_init_set(n, n_in);
    mpz_abs(n, n);                               /* allow negative n */
    if (mpz_cmp_ui(n, 1) == 0) { mpz_set_ui(order_out, 1); mpz_clear(n); return true; }

    mpz_t k, g;
    mpz_inits(k, g, NULL);
    mpz_mod(k, k_in, n);                         /* k in [0, n-1] */
    mpz_gcd(g, k, n);
    if (mpz_cmp_ui(g, 1) != 0) {
        mpz_clears(k, g, n, NULL);
        return false;                            /* not coprime */
    }
    mpz_clear(g);

    mpz_t phi;
    mpz_init(phi);
    if (!mo_eulerphi_mpz(n, phi)) {
        mpz_clears(k, phi, n, NULL);
        return false;
    }

    /* phi == 1 means the unit group is trivial; the only unit is 1 and
     * its order is 1. */
    if (mpz_cmp_ui(phi, 1) == 0) {
        mpz_set_ui(order_out, 1);
        mpz_clears(k, phi, n, NULL);
        return true;
    }

    mpz_t primes[PR_MAX_DISTINCT_PRIMES];
    size_t nprimes = 0;
    if (!pr_collect_distinct_primes(phi, primes, &nprimes)) {
        for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
        mpz_clears(k, phi, n, NULL);
        return false;
    }

    mpz_t order, q, pow;
    mpz_inits(order, q, pow, NULL);
    mpz_set(order, phi);
    for (size_t i = 0; i < nprimes; i++) {
        while (mpz_divisible_p(order, primes[i])) {
            mpz_divexact(q, order, primes[i]);
            mpz_powm(pow, k, q, n);
            if (mpz_cmp_ui(pow, 1) == 0) {
                mpz_set(order, q);
            } else {
                break;
            }
        }
    }

    mpz_set(order_out, order);
    mpz_clears(order, q, pow, NULL);
    for (size_t i = 0; i < nprimes; i++) mpz_clear(primes[i]);
    mpz_clears(k, phi, n, NULL);
    return true;
}

/* Search for the smallest m in [1, order] such that k^m mod n equals one
 * of the residues.  Pre: m_out is mpz_init'd by caller; residues are all
 * already reduced modulo n.  Returns true on success.
 *
 * Cap: bails out if order doesn't fit in unsigned long, or exceeds the
 * MO_SEARCH_CAP iteration budget — for huge groups the enumeration would
 * not terminate in a useful timeframe.  Callers see this as "no match". */
#define MO_SEARCH_CAP 100000000UL
static bool mo_search_residues(const mpz_t k_in, const mpz_t n,
                               const mpz_t order,
                               const mpz_t* residues, size_t nres,
                               mpz_t m_out) {
    if (nres == 0) return false;
    if (!mpz_fits_ulong_p(order)) return false;
    unsigned long d = mpz_get_ui(order);
    if (d > MO_SEARCH_CAP) return false;

    mpz_t cur, k;
    mpz_inits(cur, k, NULL);
    mpz_mod(k, k_in, n);
    mpz_set_ui(cur, 1);
    for (unsigned long m = 1; m <= d; m++) {
        mpz_mul(cur, cur, k);
        mpz_mod(cur, cur, n);
        for (size_t i = 0; i < nres; i++) {
            if (mpz_cmp(cur, residues[i]) == 0) {
                mpz_set_ui(m_out, m);
                mpz_clears(cur, k, NULL);
                return true;
            }
        }
    }
    mpz_clears(cur, k, NULL);
    return false;
}

static Expr* mo_emit_argt(size_t argc) {
    fprintf(stderr,
            "MultiplicativeOrder::argt: MultiplicativeOrder called with %zu "
            "argument%s; 2 or 3 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_multiplicativeorder(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return mo_emit_argt(argc);

    Expr* k_expr = res->data.function.args[0];
    Expr* n_expr = res->data.function.args[1];

    /* Integer-only contract: non-integer numerics (Real, Complex, Rational)
     * and symbolic args flow through unevaluated with no diagnostic, matching
     * MultiplicativeOrder[10., 21] -> MultiplicativeOrder[10., 21]. */
    if (!expr_is_integer_like(k_expr) || !expr_is_integer_like(n_expr)) return NULL;

    mpz_t k, n;
    expr_to_mpz(k_expr, k);
    expr_to_mpz(n_expr, n);
    if (mpz_sgn(n) == 0) { mpz_clears(k, n, NULL); return NULL; }

    mpz_t order;
    mpz_init(order);
    if (!mo_order_mpz(k, n, order)) {
        mpz_clears(k, n, order, NULL);
        return NULL;                             /* not coprime / failure */
    }

    if (argc == 2) {
        Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(order));
        mpz_clears(k, n, order, NULL);
        return out;
    }

    /* 3-arg form: third argument must be a List of integer residues.  A
     * non-List or any non-integer element leaves the call unevaluated. */
    Expr* list = res->data.function.args[2];
    if (list->type != EXPR_FUNCTION ||
        list->data.function.head->type != EXPR_SYMBOL ||
        list->data.function.head->data.symbol != SYM_List) {
        mpz_clears(k, n, order, NULL);
        return NULL;
    }
    size_t lcount = list->data.function.arg_count;
    if (lcount == 0) {
        mpz_clears(k, n, order, NULL);
        return NULL;                             /* empty target list */
    }

    mpz_t abs_n;
    mpz_init_set(abs_n, n);
    mpz_abs(abs_n, abs_n);

    mpz_t* residues = (mpz_t*)malloc(sizeof(mpz_t) * lcount);
    if (!residues) {
        mpz_clears(k, n, order, abs_n, NULL);
        return NULL;
    }
    size_t nres = 0;
    bool list_ok = true;
    for (size_t i = 0; i < lcount; i++) {
        Expr* r = list->data.function.args[i];
        if (!expr_is_integer_like(r)) { list_ok = false; break; }
        mpz_init(residues[nres]);
        expr_to_mpz(r, residues[nres]);
        mpz_mod(residues[nres], residues[nres], abs_n);
        nres++;
    }
    if (!list_ok) {
        for (size_t j = 0; j < nres; j++) mpz_clear(residues[j]);
        free(residues);
        mpz_clears(k, n, order, abs_n, NULL);
        return NULL;
    }

    mpz_t m;
    mpz_init(m);
    bool found = mo_search_residues(k, abs_n, order,
                                    (const mpz_t*)residues, nres, m);
    Expr* out = NULL;
    if (found) out = expr_bigint_normalize(expr_new_bigint_from_mpz(m));

    for (size_t i = 0; i < nres; i++) mpz_clear(residues[i]);
    free(residues);
    mpz_clears(k, n, order, abs_n, m, NULL);
    return out;
}

/* Registration, attributes, and docstrings for the number-theory builtins.
 * Called from core_init(); replaces the inline block that previously lived
 * in core.c. */
/* Emit the wrong-argument-count diagnostic for Divisible and return NULL so
 * the call is left unevaluated.  Too few arguments use Mathematica's `argm`
 * tag; too many use `argt`. */
static Expr* divisible_emit_argcount(size_t argc) {
    if (argc < 2) {
        fprintf(stderr,
                "Divisible::argm: Divisible called with %zu argument%s; "
                "2 or more arguments are expected.\n",
                argc, argc == 1 ? "" : "s");
    } else {
        fprintf(stderr,
                "Divisible::argt: Divisible called with %zu arguments; "
                "2 arguments are expected.\n",
                argc);
    }
    return NULL;
}

/* True when e is an exact Gaussian integer Complex[a, b] with both real and
 * imaginary parts integer-like (machine int or BigInt). */
static bool is_gaussian_integer(Expr* e) {
    Expr *re, *im;
    if (!is_complex(e, &re, &im)) return false;
    return expr_is_integer_like(re) && expr_is_integer_like(im);
}

/* True when e is a concrete numeric quantity (exact number, recognised
 * constant such as Pi or E, or a NumericFunction applied to numeric
 * quantities).  Mirrors the file-static helper of the same name behind
 * NumericQ in core.c — inlined here, as in rationalize.c, because that helper
 * is not exported.  Used by Divisible to distinguish a manifestly
 * non-divisible numeric pair (returns False) from symbolic arguments (left
 * unevaluated). */
static bool divisible_is_numeric_quantity(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_SYMBOL) {
        const char* name = e->data.symbol;
        return name == SYM_Pi || name == SYM_E || name == SYM_I ||
               name == SYM_Infinity || name == SYM_ComplexInfinity ||
               name == SYM_EulerGamma || name == SYM_GoldenRatio ||
               name == SYM_Catalan || name == SYM_Degree ||
               name == SYM_GoldenAngle || name == SYM_Glaisher ||
               name == SYM_Khinchin;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* head_name = e->data.function.head->data.symbol;
        if (head_name == SYM_Complex || head_name == SYM_Rational) return true;
        SymbolDef* def = symtab_get_def(head_name);
        if (def && (def->attributes & ATTR_NUMERICFUNCTION)) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (!divisible_is_numeric_quantity(e->data.function.args[i])) return false;
            }
            return true;
        }
    }
    return false;
}

/* Divisible[n, m] yields True when n is an integer multiple of m and False
 * for numeric quantities that manifestly are not.  Over the integers this is
 * Mod[n, m] == 0, computed directly with GMP's mpz_divisible_p (handling
 * BigInts and the m == 0 edge: divisible by 0 iff n == 0).  For everything
 * else — Gaussian integers, rationals, and exact numeric quantities such as
 * 2 Pi or Sqrt[6] — it forms the quotient n/m, evaluates it exactly, and
 * reports True iff that quotient collapses to an integer or Gaussian integer.
 * Symbolic, non-numeric arguments leave the call unevaluated.  Divisible is
 * Listable, so list arguments are threaded by the evaluator before reaching
 * this builtin. */
Expr* builtin_divisible(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return divisible_emit_argcount(argc);

    Expr* n = res->data.function.args[0];
    Expr* m = res->data.function.args[1];

    /* Fast path: both integer-like.  mpz_divisible_p is exact at arbitrary
     * precision, so this also covers the large-integer cases. */
    if (expr_is_integer_like(n) && expr_is_integer_like(m)) {
        mpz_t a, b;
        mpz_init(a);
        mpz_init(b);
        expr_to_mpz(n, a);
        expr_to_mpz(m, b);
        bool divisible = mpz_divisible_p(a, b) != 0;
        mpz_clears(a, b, NULL);
        return expr_new_symbol(divisible ? SYM_True : SYM_False);
    }

    /* General path: build and evaluate the quotient n * m^-1. */
    Expr* inv_args[2] = { expr_copy(m), expr_new_integer(-1) };
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), inv_args, 2);
    Expr* q_args[2] = { expr_copy(n), inv };
    Expr* q = expr_new_function(expr_new_symbol(SYM_Times), q_args, 2);
    Expr* quotient = eval_and_free(q);

    bool divisible = expr_is_integer_like(quotient) || is_gaussian_integer(quotient);
    expr_free(quotient);
    if (divisible) return expr_new_symbol(SYM_True);

    /* Not an integer.  When both arguments are concrete numeric quantities
     * they are manifestly not divisible (False); otherwise (symbolic args)
     * leave Divisible unevaluated so user rules / pattern matching apply. */
    if (divisible_is_numeric_quantity(n) && divisible_is_numeric_quantity(m)) {
        return expr_new_symbol(SYM_False);
    }
    return NULL;
}

/* ---- CoprimeQ ----------------------------------------------------------- */

/* Extract a Gaussian-integer expression into two initialized mpz_t parts.
 * Accepts a rational integer (machine int or BigInt; imaginary part 0) or an
 * exact Gaussian integer Complex[a, b] with integer-like parts.  Returns true
 * on success (re/im are set); false for anything else (rationals, reals,
 * symbols, ...).  re and im must already be mpz_init'd by the caller. */
static bool coprimeq_to_gaussian(Expr* e, mpz_t re, mpz_t im) {
    if (expr_is_integer_like(e)) {
        expr_to_mpz(e, re);
        mpz_set_ui(im, 0);
        return true;
    }
    Expr *r, *i;
    if (is_complex(e, &r, &i) && expr_is_integer_like(r) && expr_is_integer_like(i)) {
        expr_to_mpz(r, re);
        expr_to_mpz(i, im);
        return true;
    }
    return false;
}

/* Round-to-nearest integer division for the Gaussian Euclidean step:
 * q = round(num / den) with den > 0, computed as floor((2 num + den)/(2 den)).
 * All operands are mpz_t; q must be init'd by the caller. */
static void coprimeq_round_div(mpz_t q, const mpz_t num, const mpz_t den) {
    mpz_t tnum, tden;
    mpz_init(tnum);
    mpz_init(tden);
    mpz_mul_2exp(tnum, num, 1);     /* 2*num */
    mpz_add(tnum, tnum, den);       /* 2*num + den */
    mpz_mul_2exp(tden, den, 1);     /* 2*den */
    mpz_fdiv_q(q, tnum, tden);      /* floor division (den > 0) */
    mpz_clears(tnum, tden, NULL);
}

/* True iff two Gaussian integers a = a1 + a2 i and b = b1 + b2 i are coprime
 * in Z[i], i.e. their Gaussian GCD is a unit (norm 1).  Uses the Euclidean
 * algorithm with round-to-nearest division; the remainder norm strictly
 * decreases, guaranteeing termination.  Inputs are consumed by value (copied
 * internally), so the caller's mpz_t are left untouched. */
static bool gaussian_pair_coprime(const mpz_t a1, const mpz_t a2,
                                  const mpz_t b1, const mpz_t b2) {
    mpz_t x1, x2, y1, y2, q1, q2, num, t, nb, r1, r2;
    mpz_inits(x1, x2, y1, y2, q1, q2, num, t, nb, r1, r2, NULL);
    mpz_set(x1, a1); mpz_set(x2, a2);
    mpz_set(y1, b1); mpz_set(y2, b2);

    /* Loop while y != 0. */
    while (mpz_sgn(y1) != 0 || mpz_sgn(y2) != 0) {
        /* norm(y) = y1^2 + y2^2. */
        mpz_mul(nb, y1, y1);
        mpz_mul(t, y2, y2);
        mpz_add(nb, nb, t);

        /* q = round(x * conj(y) / norm(y)).
         * x*conj(y) = (x1 y1 + x2 y2) + (x2 y1 - x1 y2) i. */
        mpz_mul(num, x1, y1);
        mpz_mul(t, x2, y2);
        mpz_add(num, num, t);
        coprimeq_round_div(q1, num, nb);

        mpz_mul(num, x2, y1);
        mpz_mul(t, x1, y2);
        mpz_sub(num, num, t);
        coprimeq_round_div(q2, num, nb);

        /* r = x - q*y; q*y = (q1 y1 - q2 y2) + (q1 y2 + q2 y1) i. */
        mpz_mul(r1, q1, y1);
        mpz_mul(t, q2, y2);
        mpz_sub(r1, r1, t);
        mpz_sub(r1, x1, r1);

        mpz_mul(r2, q1, y2);
        mpz_mul(t, q2, y1);
        mpz_add(r2, r2, t);
        mpz_sub(r2, x2, r2);

        /* x <- y; y <- r. */
        mpz_set(x1, y1); mpz_set(x2, y2);
        mpz_set(y1, r1); mpz_set(y2, r2);
    }

    /* x now holds the GCD; coprime iff norm(x) == 1 (x is a unit). */
    mpz_mul(t, x1, x1);
    mpz_mul(nb, x2, x2);
    mpz_add(t, t, nb);
    bool coprime = (mpz_cmp_ui(t, 1) == 0);

    mpz_clears(x1, x2, y1, y2, q1, q2, num, t, nb, r1, r2, NULL);
    return coprime;
}

/* CoprimeQ[n1, n2, ...] yields True when the arguments are pairwise relatively
 * prime, and False otherwise.  Over the ordinary integers a pair is coprime
 * when GCD == 1; with GaussianIntegers -> True (or when any argument is an
 * exact Gaussian integer) coprimality is tested in Z[i] via the Gaussian
 * Euclidean algorithm.  Machine ints and BigInts are handled uniformly through
 * GMP.  As a *Q predicate CoprimeQ always returns a Boolean: anything that is
 * not a manifestly coprime integer or Gaussian integer (rationals, reals,
 * symbols, malformed options) makes the result False.  CoprimeQ[] is False and
 * CoprimeQ[n] (a single argument, no pairs) is True.  Listable threading over
 * list arguments is performed by the evaluator before this builtin runs. */
Expr* builtin_coprimeq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return expr_new_symbol(SYM_False);
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    /* Scan arguments: collect the numeric operands and pick up an optional
     * GaussianIntegers -> True | False rule (which, because CoprimeQ is
     * Orderless, may appear at any position).  Any other Rule, or a
     * GaussianIntegers option with a non-Boolean value, is a malformed option
     * and forces a False result. */
    bool gaussian = false;
    Expr** nums = NULL;
    size_t nnums = 0;
    if (argc > 0) nums = (Expr**)malloc(argc * sizeof(Expr*));
    for (size_t i = 0; i < argc; i++) {
        Expr* a = args[i];
        if (a->type == EXPR_FUNCTION &&
            a->data.function.head->type == EXPR_SYMBOL &&
            a->data.function.head->data.symbol == SYM_Rule &&
            a->data.function.arg_count == 2) {
            Expr* key = a->data.function.args[0];
            Expr* val = a->data.function.args[1];
            if (key->type == EXPR_SYMBOL && key->data.symbol == SYM_GaussianIntegers &&
                val->type == EXPR_SYMBOL &&
                (val->data.symbol == SYM_True || val->data.symbol == SYM_False)) {
                gaussian = (val->data.symbol == SYM_True);
                continue;
            }
            free(nums);
            return expr_new_symbol(SYM_False);  /* malformed / unknown option */
        }
        nums[nnums++] = a;
    }

    /* Per the documented contract: no operands -> False, one operand -> True. */
    if (nnums == 0) { free(nums); return expr_new_symbol(SYM_False); }
    if (nnums == 1) { free(nums); return expr_new_symbol(SYM_True); }

    /* If any operand is an exact Gaussian integer (non-real), switch to Z[i]. */
    if (!gaussian) {
        for (size_t i = 0; i < nnums; i++) {
            if (is_gaussian_integer(nums[i])) { gaussian = true; break; }
        }
    }

    bool result = true;

    if (gaussian) {
        /* Extract every operand as a Gaussian integer up front; any operand
         * that is not integer-like or a Gaussian integer is not manifestly
         * coprime -> False. */
        mpz_t* re = (mpz_t*)malloc(nnums * sizeof(mpz_t));
        mpz_t* im = (mpz_t*)malloc(nnums * sizeof(mpz_t));
        bool ok = true;
        size_t got = 0;
        for (size_t i = 0; i < nnums; i++) {
            mpz_init(re[i]);
            mpz_init(im[i]);
            got = i + 1;
            if (!coprimeq_to_gaussian(nums[i], re[i], im[i])) { ok = false; break; }
        }
        if (!ok) {
            result = false;
        } else {
            for (size_t i = 0; i < nnums && result; i++) {
                for (size_t j = i + 1; j < nnums && result; j++) {
                    if (!gaussian_pair_coprime(re[i], im[i], re[j], im[j])) result = false;
                }
            }
        }
        for (size_t i = 0; i < got; i++) mpz_clears(re[i], im[i], NULL);
        free(re);
        free(im);
    } else {
        /* Ordinary integers.  Every operand must be integer-like. */
        for (size_t i = 0; i < nnums; i++) {
            if (!expr_is_integer_like(nums[i])) { result = false; break; }
        }
        if (result) {
            mpz_t* v = (mpz_t*)malloc(nnums * sizeof(mpz_t));
            for (size_t i = 0; i < nnums; i++) {
                mpz_init(v[i]);
                expr_to_mpz(nums[i], v[i]);
            }
            mpz_t g;
            mpz_init(g);
            for (size_t i = 0; i < nnums && result; i++) {
                for (size_t j = i + 1; j < nnums && result; j++) {
                    mpz_gcd(g, v[i], v[j]);
                    if (mpz_cmp_ui(g, 1) != 0) result = false;
                }
            }
            mpz_clear(g);
            for (size_t i = 0; i < nnums; i++) mpz_clear(v[i]);
            free(v);
        }
    }

    free(nums);
    return expr_new_symbol(result ? SYM_True : SYM_False);
}

/* ---- Divisors ----------------------------------------------------------- */

/* qsort comparator for an array of mpz_t (ascending). */
static int df_mpz_cmp_qsort(const void* a, const void* b) {
    return mpz_cmp(*(const mpz_t*)a, *(const mpz_t*)b);
}

/* A Gaussian integer re + im i, used so divisor pairs can be qsort-ed as a
 * unit.  Relocating an mpz_t by memcpy (as qsort does) is safe: only the
 * struct moves, the heap-allocated limbs it points at do not. */
typedef struct { mpz_t re; mpz_t im; } df_gpair;

/* Sort Gaussian divisors by (Re, Im) ascending — matches Mathematica's
 * canonical ordering of the resulting list. */
static int df_gpair_cmp(const void* x, const void* y) {
    const df_gpair* p = (const df_gpair*)x;
    const df_gpair* q = (const df_gpair*)y;
    int c = mpz_cmp(p->re, q->re);
    if (c != 0) return c;
    return mpz_cmp(p->im, q->im);
}

/* Factor a positive integer m into distinct primes and exponents by delegating
 * to internal_factorinteger (FactorInteger).  On success returns true and sets
 * *out_primes (an array of *out_n init'd mpz_t) and *out_exps (matching
 * exponents); the caller frees both arrays and clears each prime.  For m == 1
 * the factorisation is empty (n == 0, arrays NULL).  Returns false on any
 * structural problem. */
static bool df_factor_mpz(const mpz_t m, mpz_t** out_primes,
                          unsigned long** out_exps, size_t* out_n) {
    /* internal_factorinteger takes ownership of the args array contents, so the
     * constructed integer must NOT be freed here. */
    Expr* args[1] = { expr_bigint_normalize(expr_new_bigint_from_mpz(m)) };
    Expr* fact = internal_factorinteger(args, 1);
    if (!fact || fact->type != EXPR_FUNCTION ||
        fact->data.function.head->type != EXPR_SYMBOL ||
        fact->data.function.head->data.symbol != SYM_List) {
        if (fact) expr_free(fact);
        return false;
    }

    size_t n = fact->data.function.arg_count;
    mpz_t* primes = (n > 0) ? (mpz_t*)malloc(n * sizeof(mpz_t)) : NULL;
    unsigned long* exps = (n > 0) ? (unsigned long*)malloc(n * sizeof(unsigned long)) : NULL;
    bool ok = true;
    size_t got = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* pair = fact->data.function.args[i];
        if (pair->type != EXPR_FUNCTION || pair->data.function.arg_count != 2 ||
            !expr_is_integer_like(pair->data.function.args[0]) ||
            pair->data.function.args[1]->type != EXPR_INTEGER ||
            pair->data.function.args[1]->data.integer <= 0) {
            ok = false; break;
        }
        expr_to_mpz(pair->data.function.args[0], primes[i]);  /* inits primes[i] */
        exps[i] = (unsigned long)pair->data.function.args[1]->data.integer;
        got = i + 1;
    }
    expr_free(fact);

    if (!ok) {
        for (size_t i = 0; i < got; i++) mpz_clear(primes[i]);
        free(primes);
        free(exps);
        return false;
    }
    *out_primes = primes;
    *out_exps = exps;
    *out_n = n;
    return true;
}

/* Build the ascending list of positive divisors of nabs (>= 2 in practice).
 * Enumerates the divisor lattice with a mixed-radix exponent odometer over the
 * prime factorisation.  Returns a List Expr, or NULL if factoring fails. */
static Expr* divisors_ordinary(const mpz_t nabs) {
    mpz_t* primes;
    unsigned long* exps;
    size_t np;
    if (!df_factor_mpz(nabs, &primes, &exps, &np)) return NULL;

    /* Number of divisors = prod(exps[i] + 1).  Bail (leave unevaluated) if this
     * product overflows size_t — e.g. Divisors[100!] has ~10^28 divisors, which
     * is intractable to materialise as a list. */
    size_t ndiv = 1;
    bool overflow = false;
    for (size_t i = 0; i < np; i++) {
        size_t f = exps[i] + 1;
        if (ndiv > SIZE_MAX / f) { overflow = true; break; }
        ndiv *= f;
    }
    if (overflow) {
        for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
        free(primes);
        free(exps);
        return NULL;
    }

    mpz_t* divs = (mpz_t*)malloc(ndiv * sizeof(mpz_t));
    if (!divs) {
        for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
        free(primes);
        free(exps);
        return NULL;
    }
    unsigned long* cur = (unsigned long*)calloc(np, sizeof(unsigned long));
    mpz_t acc, pw;
    mpz_inits(acc, pw, NULL);
    for (size_t d = 0; d < ndiv; d++) {
        mpz_set_ui(acc, 1);
        for (size_t i = 0; i < np; i++) {
            mpz_pow_ui(pw, primes[i], cur[i]);
            mpz_mul(acc, acc, pw);
        }
        mpz_init_set(divs[d], acc);
        /* advance mixed-radix odometer: digit i ranges 0..exps[i] */
        for (size_t i = 0; i < np; i++) {
            if (++cur[i] <= exps[i]) break;
            cur[i] = 0;
        }
    }
    mpz_clears(acc, pw, NULL);
    free(cur);

    qsort(divs, ndiv, sizeof(mpz_t), df_mpz_cmp_qsort);

    Expr** items = (Expr**)malloc(ndiv * sizeof(Expr*));
    for (size_t d = 0; d < ndiv; d++)
        items[d] = expr_bigint_normalize(expr_new_bigint_from_mpz(divs[d]));
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, ndiv);
    free(items);

    for (size_t d = 0; d < ndiv; d++) mpz_clear(divs[d]);
    free(divs);
    for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
    free(primes);
    free(exps);
    return out;
}

/* Solve p = c^2 + d^2 (c, d > 0) for a prime p ≡ 1 (mod 4) via Cornacchia's
 * algorithm.  Finds a square root of -1 modulo p (t^((p-1)/4) for a quadratic
 * non-residue t), then Euclidean-reduces (p, root) until the remainder drops
 * below sqrt(p).  Sets c, d (caller pre-inits) and returns true on success. */
static bool df_two_squares(const mpz_t p, mpz_t c, mpz_t d) {
    mpz_t t, x, e, a, b, r, s, sq;
    mpz_inits(t, x, e, a, b, r, s, sq, NULL);

    bool found = false;
    for (unsigned long tv = 2; tv < 1000000UL; tv++) {
        mpz_set_ui(t, tv);
        if (mpz_legendre(t, p) == -1) { found = true; break; }
    }
    bool ok = false;
    if (found) {
        mpz_sub_ui(e, p, 1);
        mpz_fdiv_q_ui(e, e, 4);
        mpz_powm(x, t, e, p);           /* x^2 ≡ -1 (mod p) */

        mpz_set(a, p);
        mpz_set(b, x);
        mpz_mul(sq, b, b);
        while (mpz_cmp(sq, p) >= 0) {    /* reduce until b^2 < p */
            mpz_mod(r, a, b);
            mpz_set(a, b);
            mpz_set(b, r);
            mpz_mul(sq, b, b);
        }
        mpz_set(c, b);
        mpz_mul(sq, b, b);
        mpz_sub(s, p, sq);              /* s = p - c^2 */
        mpz_sqrt(d, s);
        mpz_mul(sq, d, d);
        ok = (mpz_cmp(sq, s) == 0);     /* confirm s is a perfect square */
    }

    mpz_clears(t, x, e, a, b, r, s, sq, NULL);
    return ok;
}

/* If the Gaussian integer g = u + v i divides z = a + b i exactly, overwrite
 * (a, b) with the quotient and return true; otherwise leave (a, b) unchanged
 * and return false.  z/g = z conj(g) / N(g). */
static bool df_gauss_divide(mpz_t a, mpz_t b, const mpz_t u, const mpz_t v) {
    mpz_t nrm, re, im, t;
    mpz_inits(nrm, re, im, t, NULL);
    mpz_mul(nrm, u, u); mpz_mul(t, v, v); mpz_add(nrm, nrm, t);   /* N(g) */
    mpz_mul(re, a, u); mpz_mul(t, b, v); mpz_add(re, re, t);      /* a u + b v */
    mpz_mul(im, b, u); mpz_mul(t, a, v); mpz_sub(im, im, t);      /* b u - a v */
    bool ok = mpz_divisible_p(re, nrm) && mpz_divisible_p(im, nrm);
    if (ok) {
        mpz_divexact(a, re, nrm);
        mpz_divexact(b, im, nrm);
    }
    mpz_clears(nrm, re, im, t, NULL);
    return ok;
}

/* Rotate (re, im) by multiples of i into the canonical first-quadrant
 * associate: the unique unit-multiple with Re > 0 and Im >= 0. */
static void df_normalize_quadrant(mpz_t re, mpz_t im) {
    mpz_t t;
    mpz_init(t);
    for (int k = 0; k < 4; k++) {
        if (mpz_sgn(re) > 0 && mpz_sgn(im) >= 0) break;
        mpz_neg(t, im);    /* multiply by i: (re, im) -> (-im, re) */
        mpz_set(im, re);
        mpz_set(re, t);
    }
    mpz_clear(t);
}

/* Extract a Gaussian-integer expression into two mpz_t parts a, b (both are
 * INITIALIZED here — the caller must NOT pre-init them).  Accepts a rational
 * integer (machine int / BigInt; im 0) or Complex[a, b] with integer-like
 * parts.  Returns false (leaving a, b uninitialised) for anything else. */
static bool df_to_gaussian(Expr* e, mpz_t a, mpz_t b) {
    if (expr_is_integer_like(e)) {
        expr_to_mpz(e, a);          /* inits a */
        mpz_init_set_ui(b, 0);
        return true;
    }
    Expr *r, *i;
    if (is_complex(e, &r, &i) && expr_is_integer_like(r) && expr_is_integer_like(i)) {
        expr_to_mpz(r, a);          /* inits a */
        expr_to_mpz(i, b);          /* inits b */
        return true;
    }
    return false;
}

/* Build the list of Gaussian divisors of z = a_in + b_in i, one representative
 * per associate class (first quadrant), sorted by (Re, Im).  Factors N(z) over
 * Z, lifts each rational prime to its Gaussian prime(s), enumerates the divisor
 * lattice, and normalises each product to its first-quadrant associate.
 * Returns a List Expr, or NULL on z == 0 / factoring failure. */
static Expr* divisors_gaussian(const mpz_t a_in, const mpz_t b_in) {
    mpz_t N, t;
    mpz_inits(N, t, NULL);
    mpz_mul(N, a_in, a_in);
    mpz_mul(t, b_in, b_in);
    mpz_add(N, N, t);
    mpz_clear(t);

    if (mpz_sgn(N) == 0) { mpz_clear(N); return NULL; }   /* z = 0 */

    if (mpz_cmp_ui(N, 1) == 0) {                           /* unit: only 1 */
        mpz_clear(N);
        Expr* items[1] = { expr_new_integer(1) };
        return expr_new_function(expr_new_symbol(SYM_List), items, 1);
    }

    mpz_t* primes;
    unsigned long* exps;
    size_t np;
    bool fok = df_factor_mpz(N, &primes, &exps, &np);
    mpz_clear(N);
    if (!fok) return NULL;

    /* Gaussian prime factors of z: (gu[i] + gv[i] i)^ge[i]. */
    size_t cap = 2 * np + 1;
    mpz_t* gu = (mpz_t*)malloc(cap * sizeof(mpz_t));
    mpz_t* gv = (mpz_t*)malloc(cap * sizeof(mpz_t));
    unsigned long* ge = (unsigned long*)malloc(cap * sizeof(unsigned long));
    size_t gn = 0;
    bool ok = true;

    for (size_t i = 0; i < np && ok; i++) {
        if (mpz_cmp_ui(primes[i], 2) == 0) {
            /* 2 ramifies: (1 + i) with exponent v_2(N). */
            mpz_init_set_ui(gu[gn], 1);
            mpz_init_set_ui(gv[gn], 1);
            ge[gn] = exps[i];
            gn++;
        } else if (mpz_fdiv_ui(primes[i], 4) == 3) {
            /* p ≡ 3 (mod 4) is inert; v_p(N) is even, exponent in z is half. */
            mpz_init_set(gu[gn], primes[i]);
            mpz_init_set_ui(gv[gn], 0);
            ge[gn] = exps[i] / 2;
            gn++;
        } else {
            /* p ≡ 1 (mod 4) splits as (c + d i)(c - d i). */
            mpz_t c, d;
            mpz_inits(c, d, NULL);
            if (!df_two_squares(primes[i], c, d)) {
                mpz_clears(c, d, NULL);
                ok = false;
                break;
            }
            /* Exponent of c + d i in z by repeated exact division. */
            mpz_t wa, wb;
            mpz_init_set(wa, a_in);
            mpz_init_set(wb, b_in);
            unsigned long epi = 0;
            while (df_gauss_divide(wa, wb, c, d)) epi++;
            mpz_clears(wa, wb, NULL);
            unsigned long econj = exps[i] - epi;
            if (epi > 0) {
                mpz_init_set(gu[gn], c);
                mpz_init_set(gv[gn], d);
                ge[gn] = epi;
                gn++;
            }
            if (econj > 0) {
                mpz_init_set(gu[gn], c);
                mpz_init(gv[gn]);
                mpz_neg(gv[gn], d);
                ge[gn] = econj;
                gn++;
            }
            mpz_clears(c, d, NULL);
        }
    }

    for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
    free(primes);
    free(exps);

    if (!ok) {
        for (size_t i = 0; i < gn; i++) mpz_clears(gu[i], gv[i], NULL);
        free(gu); free(gv); free(ge);
        return NULL;
    }

    size_t ndiv = 1;
    bool ovf = false;
    for (size_t i = 0; i < gn; i++) {
        size_t f = ge[i] + 1;
        if (ndiv > SIZE_MAX / f) { ovf = true; break; }
        ndiv *= f;
    }
    df_gpair* divs = ovf ? NULL : (df_gpair*)malloc(ndiv * sizeof(df_gpair));
    if (!divs) {
        for (size_t i = 0; i < gn; i++) mpz_clears(gu[i], gv[i], NULL);
        free(gu); free(gv); free(ge);
        return NULL;
    }
    unsigned long* cur = (unsigned long*)calloc(gn, sizeof(unsigned long));
    mpz_t re, im, pu, pv, nu, nv, tt;
    mpz_inits(re, im, pu, pv, nu, nv, tt, NULL);
    for (size_t d = 0; d < ndiv; d++) {
        mpz_set_ui(re, 1);
        mpz_set_ui(im, 0);
        for (size_t i = 0; i < gn; i++) {
            /* (pu + pv i) = (gu[i] + gv[i] i)^cur[i] */
            mpz_set_ui(pu, 1);
            mpz_set_ui(pv, 0);
            for (unsigned long k = 0; k < cur[i]; k++) {
                mpz_mul(nu, pu, gu[i]); mpz_mul(tt, pv, gv[i]); mpz_sub(nu, nu, tt);
                mpz_mul(nv, pu, gv[i]); mpz_mul(tt, pv, gu[i]); mpz_add(nv, nv, tt);
                mpz_set(pu, nu); mpz_set(pv, nv);
            }
            /* (re + im i) *= (pu + pv i) */
            mpz_mul(nu, re, pu); mpz_mul(tt, im, pv); mpz_sub(nu, nu, tt);
            mpz_mul(nv, re, pv); mpz_mul(tt, im, pu); mpz_add(nv, nv, tt);
            mpz_set(re, nu); mpz_set(im, nv);
        }
        df_normalize_quadrant(re, im);
        mpz_init_set(divs[d].re, re);
        mpz_init_set(divs[d].im, im);
        for (size_t i = 0; i < gn; i++) {
            if (++cur[i] <= ge[i]) break;
            cur[i] = 0;
        }
    }
    mpz_clears(re, im, pu, pv, nu, nv, tt, NULL);
    free(cur);
    for (size_t i = 0; i < gn; i++) mpz_clears(gu[i], gv[i], NULL);
    free(gu); free(gv); free(ge);

    qsort(divs, ndiv, sizeof(df_gpair), df_gpair_cmp);

    Expr** items = (Expr**)malloc(ndiv * sizeof(Expr*));
    for (size_t d = 0; d < ndiv; d++) {
        if (mpz_sgn(divs[d].im) == 0) {
            items[d] = expr_bigint_normalize(expr_new_bigint_from_mpz(divs[d].re));
        } else {
            Expr* cargs[2] = {
                expr_bigint_normalize(expr_new_bigint_from_mpz(divs[d].re)),
                expr_bigint_normalize(expr_new_bigint_from_mpz(divs[d].im))
            };
            items[d] = expr_new_function(expr_new_symbol(SYM_Complex), cargs, 2);
        }
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, ndiv);
    free(items);
    for (size_t d = 0; d < ndiv; d++) mpz_clears(divs[d].re, divs[d].im, NULL);
    free(divs);
    return out;
}

/* Emit the wrong-argument-count diagnostic for Divisors and return NULL so the
 * call is left unevaluated. */
static Expr* divisors_emit_argx(size_t npos) {
    if (!arith_warnings_muted()) {
        fprintf(stderr,
                "Divisors::argx: Divisors called with %zu argument%s; "
                "1 argument is expected.\n",
                npos, npos == 1 ? "" : "s");
    }
    return NULL;
}

/* Divisors[n] gives the ascending list of positive integers dividing n (sign
 * ignored).  Divisors[n, GaussianIntegers -> True] (or a non-real Gaussian
 * input) returns the divisors in Z[i], one first-quadrant representative per
 * associate class, sorted by (Re, Im).  Machine ints and BigInts are handled
 * uniformly through GMP.  Divisors is Listable, so list arguments are threaded
 * by the evaluator before this builtin runs.  Non-integer, zero, or 0-argument
 * calls are left unevaluated (the last with a Divisors::argx message). */
Expr* builtin_divisors(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    /* Separate the lone positional argument from a GaussianIntegers option. */
    bool gaussian = false, gaussian_set = false;
    Expr* n = NULL;
    size_t npos = 0;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = args[i];
        if (a->type == EXPR_FUNCTION &&
            a->data.function.head->type == EXPR_SYMBOL &&
            a->data.function.head->data.symbol == SYM_Rule &&
            a->data.function.arg_count == 2) {
            Expr* key = a->data.function.args[0];
            Expr* val = a->data.function.args[1];
            if (key->type == EXPR_SYMBOL && key->data.symbol == SYM_GaussianIntegers &&
                val->type == EXPR_SYMBOL &&
                (val->data.symbol == SYM_True || val->data.symbol == SYM_False)) {
                gaussian = (val->data.symbol == SYM_True);
                gaussian_set = true;
                continue;
            }
            return NULL;  /* malformed / unknown option: leave unevaluated */
        }
        npos++;
        n = a;
    }

    if (npos != 1) return divisors_emit_argx(npos);

    /* Explicit option wins; otherwise auto-enable Z[i] for non-real input. */
    bool use_gaussian = gaussian_set ? gaussian : is_gaussian_integer(n);

    if (use_gaussian) {
        mpz_t a, b;
        if (!df_to_gaussian(n, a, b)) return NULL;
        Expr* out = divisors_gaussian(a, b);
        mpz_clears(a, b, NULL);
        return out;
    }

    /* Ordinary integer path. */
    if (!expr_is_integer_like(n)) return NULL;
    mpz_t v;
    expr_to_mpz(n, v);                 /* inits v */
    if (mpz_sgn(v) == 0) { mpz_clear(v); return NULL; }   /* Divisors[0] */
    mpz_abs(v, v);

    Expr* out;
    if (mpz_cmp_ui(v, 1) == 0) {
        Expr* items[1] = { expr_new_integer(1) };
        out = expr_new_function(expr_new_symbol(SYM_List), items, 1);
    } else {
        out = divisors_ordinary(v);
    }
    mpz_clear(v);
    return out;
}

void numbertheory_init(void) {
    symtab_add_builtin("GCD", builtin_gcd);
    symtab_add_builtin("LCM", builtin_lcm);
    symtab_add_builtin("ExtendedGCD", builtin_extendedgcd);
    symtab_add_builtin("PowerMod", builtin_powermod);
    symtab_add_builtin("PrimitiveRoot", builtin_primitiveroot);
    symtab_add_builtin("PrimitiveRootList", builtin_primitiverootlist);
    symtab_add_builtin("MultiplicativeOrder", builtin_multiplicativeorder);
    symtab_add_builtin("Factorial", builtin_factorial);
    symtab_add_builtin("Factorial2", builtin_factorial2);
    symtab_add_builtin("FactorialPower", builtin_factorialpower);
    symtab_add_builtin("Binomial", builtin_binomial);
    symtab_add_builtin("Divisible", builtin_divisible);
    symtab_add_builtin("CoprimeQ", builtin_coprimeq);
    symtab_add_builtin("Divisors", builtin_divisors);

    symtab_get_def("GCD")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY);
    symtab_get_def("LCM")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY);
    symtab_get_def("ExtendedGCD")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    symtab_get_def("PowerMod")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_get_def("PrimitiveRoot")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_get_def("PrimitiveRootList")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_get_def("MultiplicativeOrder")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Factorial")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Factorial2")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("FactorialPower")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_set_docstring("FactorialPower",
        "FactorialPower[n, k]\n\tThe falling factorial n (n - 1) (n - 2) ... (n - k + 1).\n\tFor non-negative integer k, expands to a product of k linear factors.\n\tEquivalent to n! / (n - k)! when both n and k are non-negative integers.");
    symtab_get_def("Binomial")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Divisible")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    symtab_get_def("CoprimeQ")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE | ATTR_ORDERLESS);
    symtab_get_def("Divisors")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
}


#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

int g_arith_warnings_muted = 0;

int64_t gcd(int64_t a, int64_t b) {
    a = llabs(a);
    b = llabs(b);
    while (b) {
        a %= b;
        int64_t tmp = a;
        a = b;
        b = tmp;
    }
    return a;
}

int64_t lcm(int64_t a, int64_t b) {
    if (a == 0 || b == 0) return 0;
    a = llabs(a);
    b = llabs(b);
    return (a / gcd(a, b)) * b;
}

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
        out = expr_new_function(expr_new_symbol("Rational"), args, 2);
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

Expr* make_rational(int64_t n, int64_t d) {
    if (d == 0) return NULL; // Error
    if (n == 0) return expr_new_integer(0);
    
    int64_t common = gcd(n, d);
    n /= common;
    d /= common;

    if (d < 0) {
        n = -n;
        d = -d;
    }

    if (d == 1) return expr_new_integer(n);

    Expr* args[2];
    args[0] = expr_new_integer(n);
    args[1] = expr_new_integer(d);
    return expr_new_function(expr_new_symbol("Rational"), args, 2);
}

Expr* builtin_rational(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* n_expr = res->data.function.args[0];
    Expr* d_expr = res->data.function.args[1];
    
    if (n_expr->type == EXPR_INTEGER && d_expr->type == EXPR_INTEGER) {
        int64_t n = n_expr->data.integer;
        int64_t d = d_expr->data.integer;
        if (d == 0) {
            if (!arith_warnings_muted())
                fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
            if (n == 0) {
                if (!arith_warnings_muted())
                    fprintf(stderr,
                        "Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.\n");
                return expr_new_symbol("Indeterminate");
            }
            return expr_new_symbol("ComplexInfinity");
        }
        
        Expr* r = make_rational(n, d);
        if (r && r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL && r->data.function.head->data.symbol == SYM_Rational) {
            Expr* rn = r->data.function.args[0];
            Expr* rd = r->data.function.args[1];
            if (rn->type == EXPR_INTEGER && rd->type == EXPR_INTEGER && rn->data.integer == n && rd->data.integer == d) {
                // No simplification happened
                expr_free(r);
                return NULL;
            }
        }
        return r;
    }
    return NULL;
}

bool is_rational(const Expr* e, int64_t* n, int64_t* d) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) {
        if (n) *n = e->data.integer;
        if (d) *d = 1;
        return true;
    }
    /* Function nodes built during evaluator transitions can transiently
     * carry a NULL head pointer; guard before dereferencing.  Without
     * this check, a Function with a NULL head reaching this fast-path
     * crashes the dispatcher (observed on Linux under nested-radical
     * Simplify; macOS heap layout makes the deref usually land in
     * mapped memory). */
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Rational) {
        if (e->data.function.arg_count == 2 &&
            e->data.function.args[0] &&
            e->data.function.args[1] &&
            e->data.function.args[0]->type == EXPR_INTEGER &&
            e->data.function.args[1]->type == EXPR_INTEGER) {
            if (n) *n = e->data.function.args[0]->data.integer;
            if (d) *d = e->data.function.args[1]->data.integer;
            return true;
        }
    }
    return false;
}

/* Bigint-aware version: matches Integer, BigInt, or Rational[X, Y] where
 * both components are integer-like (Integer or BigInt).  Used in numeric
 * predicates so the Times/Plus folders recognise rationals whose
 * components have overflowed int64 — without this, a BigInt times a
 * Rational[1, BigInt] is left unsimplified and can defeat exact
 * polynomial division. */
bool is_rational_like(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Rational &&
        e->data.function.arg_count == 2 &&
        expr_is_integer_like(e->data.function.args[0]) &&
        expr_is_integer_like(e->data.function.args[1])) {
        return true;
    }
    return false;
}

bool is_complex(Expr* e, Expr** re, Expr** im) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Complex &&
        e->data.function.arg_count == 2) {
        if (re) *re = e->data.function.args[0];
        if (im) *im = e->data.function.args[1];
        return true;
    }
    return false;
}

Expr* make_complex(Expr* re, Expr* im) {
    Expr* args[2] = { re, im };
    return expr_new_function(expr_new_symbol("Complex"), args, 2);
}

Expr* builtin_subtract(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];

    Expr* minus_one = expr_new_integer(-1);
    Expr* mb_args[2] = { minus_one, expr_copy(b) };
    Expr* minus_b = expr_new_function(expr_new_symbol("Times"), mb_args, 2);

    Expr* p_args[2] = { expr_copy(a), minus_b };
    return expr_new_function(expr_new_symbol("Plus"), p_args, 2);
}

Expr* builtin_complex(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* re = res->data.function.args[0];
    Expr* im = res->data.function.args[1];

    if (im->type == EXPR_INTEGER && im->data.integer == 0) {
        return expr_copy(re);
    }
    if (im->type == EXPR_REAL && im->data.real == 0.0) {
        if (re->type == EXPR_INTEGER) {
            return expr_new_real((double)re->data.integer);
        }
        return expr_copy(re);
    }

    return NULL;
}

Expr* builtin_divide(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    
    Expr* num = res->data.function.args[0];
    Expr* den = res->data.function.args[1];

    if (num->type == EXPR_REAL || den->type == EXPR_REAL) {
        double vnum = (num->type == EXPR_REAL) ? num->data.real : (num->type == EXPR_INTEGER) ? (double)num->data.integer : (num->type == EXPR_BIGINT) ? mpz_get_d(num->data.bigint) : 0.0;
        double vden = (den->type == EXPR_REAL) ? den->data.real : (den->type == EXPR_INTEGER) ? (double)den->data.integer : (den->type == EXPR_BIGINT) ? mpz_get_d(den->data.bigint) : 0.0;
        if (vden == 0.0) {
            if (!arith_warnings_muted())
                fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
            return expr_new_symbol("ComplexInfinity");
        }
        return expr_new_real(vnum / vden);
    }

    int64_t n1, d1, n2, d2;
    if (is_rational(num, &n1, &d1) && is_rational(den, &n2, &d2)) {
        if (n2 == 0) {
            /* x / 0 with rational/integer x: 0/0 -> Indeterminate (handled in
             * Times when 0 multiplies ComplexInfinity); otherwise emit the
             * Power::infy message and yield ComplexInfinity. */
            if (!arith_warnings_muted())
                fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
            if (n1 == 0) {
                if (!arith_warnings_muted())
                    fprintf(stderr,
                        "Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.\n");
                return expr_new_symbol("Indeterminate");
            }
            return expr_new_symbol("ComplexInfinity");
        }
        Expr* r = make_rational(n1 * d2, d1 * n2);
        if (r) return r;
    }

    Expr* minus_one = expr_new_integer(-1);
    Expr* p_args[2] = { expr_copy(den), minus_one };
    Expr* power = expr_new_function(expr_new_symbol("Power"), p_args, 2);
    
    Expr* t_args[2] = { expr_copy(num), power };
    return expr_new_function(expr_new_symbol("Times"), t_args, 2);
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
            if (n < 0) return expr_new_symbol("ComplexInfinity");
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
            int64_t num = 1;
            int64_t den = 1;
            
            if (n > 0) {
                for (int64_t i = n; i >= 1; i -= 2) num *= i;
                den = 1LL << ((n + 1) / 2);
            } else {
                for (int64_t i = n + 2; i <= -1; i += 2) {
                    num *= 2;
                    den *= i;
                }
            }
            
            Expr* coeff = make_rational(num, den);
            if (!coeff) coeff = expr_new_integer(0);
            
            Expr* pi_sym = expr_new_symbol("Pi");
            Expr* half = make_rational(1, 2);
            Expr* sqrt_pi = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){pi_sym, half}, 2));
            
            if (coeff->type == EXPR_INTEGER && coeff->data.integer == 1) {
                expr_free(coeff);
                return sqrt_pi;
            } else {
                return eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){coeff, sqrt_pi}, 2));
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
                factors[i] = expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(n_arg), shift }, 2);
            }
        }
        Expr* product = expr_new_function(expr_new_symbol("Times"), factors, (size_t)k);
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
            return expr_new_symbol("ComplexInfinity");
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
            return expr_new_symbol("ComplexInfinity");
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
            factors[i] = expr_new_function(expr_new_symbol("Plus"),
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
        inv = expr_new_function(expr_new_symbol("Rational"),
            (Expr*[]){ one, fact_expr }, 2);
    }
    factors[k] = inv;

    Expr* product = expr_new_function(expr_new_symbol("Times"),
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
        Expr* diff = expr_new_function(expr_new_symbol("Subtract"),
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

bool is_infinity_sym(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_Infinity;
}

bool is_complex_infinity_sym(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_ComplexInfinity;
}

bool is_indeterminate_sym(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_Indeterminate;
}

int expr_numeric_sign(Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer > 0) return 1;
        if (e->data.integer < 0) return -1;
        return 0;
    }
    if (e->type == EXPR_REAL) {
        if (e->data.real > 0.0) return 1;
        if (e->data.real < 0.0) return -1;
        return 0;
    }
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint);
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        // d is conventionally positive in Mathilda Rational[n, d].
        if (n > 0) return (d > 0) ? 1 : -1;
        if (n < 0) return (d > 0) ? -1 : 1;
        return 0;
    }
    return 0;
}

bool expr_is_superficially_negative(Expr* e) {
    if (!e) return false;
    int s = expr_numeric_sign(e);
    if (s < 0) return true;
    if (s > 0) return false;
    /* s == 0: either a zero numeric (not negative) or a non-numeric --
     * fall through to Complex and Times shape checks. */
    Expr* re; Expr* im;
    if (is_complex(e, &re, &im)) {
        int rs = expr_numeric_sign(re);
        if (rs < 0) return true;
        if (rs > 0) return false;
        /* Real part is zero: use the imaginary part's sign. Catches the
         * "pure imaginary with negative coefficient" case, e.g. -2 I =
         * Complex[0, -2]. */
        return expr_numeric_sign(im) < 0;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Times &&
        e->data.function.arg_count > 0) {
        /* Leading factor carries the syntactic sign, per Times canonical
         * ordering (numerics sort first). */
        return expr_is_superficially_negative(e->data.function.args[0]);
    }
    return false;
}

bool is_neg_infinity_form(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    if (e->data.function.head->data.symbol != SYM_Times) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (!is_infinity_sym(b)) return false;
    return expr_numeric_sign(a) < 0;
}

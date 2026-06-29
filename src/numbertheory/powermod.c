/* powermod.c -- PowerMod[] and its modular-root helpers.
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

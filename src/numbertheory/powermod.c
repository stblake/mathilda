/* powermod.c -- PowerMod[] and its modular-root helpers.
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "ndarray.h"      /* is_ndarray, ndarray_delist_and_reeval */
#include "ndinteger.h"   /* NDArray buffer fast path */
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "common.h"      /* head_is */
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "message.h"   /* mth_message_gated: Quiet/Check funnel */
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
 * r-th roots, so we build one that is COMPLETE -- it either finds a root,
 * proves there is none, or (only for a root degree with a prime factor
 * too large for the discrete-log step) reports the case unsupported:
 *
 *   1. Factor m (a prime is taken as is) and solve per prime power p^e.
 *   2. For x^r ≡ c (mod p^e): c ≡ 0 gives x = 0; otherwise write
 *      c = p^v u with u a unit.  A root exists only when r | v, and then
 *      x = p^(v/r) y with y^r ≡ u (mod p^(e-v)).
 *   3. The unit root works in the unit group, which is cyclic of order
 *      p^(k-1)(p-1) for odd p (a square root there keeps the established
 *      Tonelli-Shanks + Hensel lift, valid since p ∤ 2; a prime power up to
 *      MODROOT_BRUTE_LIMIT is searched exhaustively for its smallest root).
 *      For p = 2 (k >= 3) the group is
 *      {±1} x <5>: an odd r is inverted modulo the exponent 2^(k-2); an
 *      even r needs u ≡ 1 (mod 4) and is solved in the cyclic <5>.
 *   4. In a cyclic group of order N, x^r = u is solvable iff
 *      u^(N/g) = 1 with g = gcd(r, N).  Since gcd(r/g, N/g) = 1,
 *      y = u^((r/g)^-1 mod N/g) is a g-th power with y^(r/g) = u, and the
 *      g-th root of y is taken one prime q | g at a time with the
 *      Adleman-Manders-Miller (generalised Tonelli-Shanks) step.  Any q-th
 *      root of a g-th power is again a (g/q)-th power (mu_q lies in the
 *      subgroup of (g/q)-th powers because g | N), so no branch fails.
 *   5. Combine the per-prime-power roots with CRT, and verify x^r ≡ c.
 *
 * Earlier versions Hensel-lifted a root mod p, which is invalid when p | r
 * (square roots mod 2^k) or when the root is ≡ 0 mod p, and used brute force
 * for r > 2 with gcd(r, p - 1) > 1; each of those reported "no solution"
 * for solvable equations (PowerMod#1-#3 of the 2026-11 audit).
 */

#define MODROOT_BRUTE_LIMIT 1000000UL
/* Largest prime root-degree factor q whose discrete log (in a group of order
 * q) the AMM step searches by enumeration. */
#define MODROOT_DLOG_LIMIT 1048576UL

/* Outcome of a modular-root search. */
enum { MR_UNSUPPORTED = -1, MR_NONE = 0, MR_FOUND = 1 };

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

/* q-th root (q prime, q | N) of y in a cyclic group of order N inside
 * (Z/mod)^*, for y a q-th power (y^(N/q) = 1).  z_hint, when given, is an
 * element of the group that is not a q-th power; otherwise one is searched
 * for.  Adleman-Manders-Miller: N = q^t s with q ∤ s; x = y^(q^-1 mod s)
 * leaves an error e = x^q / y in the q-Sylow subgroup, which is cancelled
 * one q-adic digit at a time with powers of c = z^s (order q^t).  x must not
 * alias y. */
static int amm_root(mpz_t x, const mpz_t y, unsigned long q, const mpz_t N,
                    const mpz_t mod, const mpz_t z_hint) {
    mpz_t s, Nq, z, c, zeta, e, f, w, tmp, qq, alpha;
    mpz_inits(s, Nq, z, c, zeta, e, f, w, tmp, qq, alpha, NULL);
    int rc = MR_FOUND;
    unsigned long t = 0;
    mpz_set(s, N);
    while (mpz_divisible_ui_p(s, q)) { mpz_divexact_ui(s, s, q); t++; }
    mpz_divexact_ui(Nq, N, q);
    mpz_set_ui(qq, q);

    if (z_hint) {
        mpz_set(z, z_hint);
    } else {
        bool found = false;
        for (unsigned long k = 2; k < 100000UL && !found; k++) {
            mpz_set_ui(z, k);
            if (mpz_cmp(z, mod) >= 0) break;
            mpz_gcd(tmp, z, mod);
            if (mpz_cmp_ui(tmp, 1) != 0) continue;
            mpz_powm(tmp, z, Nq, mod);
            found = mpz_cmp_ui(tmp, 1) != 0;
        }
        if (!found) { rc = MR_UNSUPPORTED; goto done; }
    }
    mpz_powm(c, z, s, mod);                          /* order exactly q^t */

    if (mpz_cmp_ui(s, 1) == 0) mpz_set_ui(alpha, 0);
    else mpz_invert(alpha, qq, s);
    mpz_powm(x, y, alpha, mod);
    mpz_powm(e, x, qq, mod);                         /* e = x^q / y */
    if (!mpz_invert(tmp, y, mod)) { rc = MR_UNSUPPORTED; goto done; }
    mpz_mul(e, e, tmp);
    mpz_mod(e, e, mod);

    mpz_set(zeta, c);                                /* zeta = c^(q^(t-1)) */
    for (unsigned long k = 1; k < t; k++) mpz_powm(zeta, zeta, qq, mod);

    while (mpz_cmp_ui(e, 1) != 0) {
        /* Smallest i >= 1 with e^(q^i) = 1; f = e^(q^(i-1)) has order q. */
        unsigned long i = 1;
        mpz_set(f, e);
        for (;;) {
            mpz_powm(tmp, f, qq, mod);
            if (mpz_cmp_ui(tmp, 1) == 0) break;
            mpz_set(f, tmp);
            if (++i >= t) { rc = MR_UNSUPPORTED; goto done; }  /* y not a q-th power */
        }
        if (q > MODROOT_DLOG_LIMIT) { rc = MR_UNSUPPORTED; goto done; }
        unsigned long j = 0;                         /* zeta^j = f */
        mpz_set_ui(tmp, 1);
        for (unsigned long k = 1; k < q; k++) {
            mpz_mul(tmp, tmp, zeta);
            mpz_mod(tmp, tmp, mod);
            if (mpz_cmp(tmp, f) == 0) { j = k; break; }
        }
        if (j == 0) { rc = MR_UNSUPPORTED; goto done; }
        /* w = c^(j q^(t-1-i)); x <- x / w, e <- e / w^q. */
        mpz_set(w, c);
        for (unsigned long k = 0; k + i + 1 < t; k++) mpz_powm(w, w, qq, mod);
        mpz_powm_ui(w, w, j, mod);
        if (!mpz_invert(w, w, mod)) { rc = MR_UNSUPPORTED; goto done; }
        mpz_mul(x, x, w);
        mpz_mod(x, x, mod);
        mpz_powm(tmp, w, qq, mod);
        mpz_mul(e, e, tmp);
        mpz_mod(e, e, mod);
    }
done:
    mpz_clears(s, Nq, z, c, zeta, e, f, w, tmp, qq, alpha, NULL);
    return rc;
}

/* x with x^r = u in a cyclic group of order N inside (Z/mod)^* (u in the
 * group).  MR_NONE when u is not an r-th power there.  See step 4 above. */
static int cyc_root(mpz_t x, const mpz_t u, const mpz_t r, const mpz_t N,
                    const mpz_t mod, const mpz_t z_hint) {
    mpz_t g, Ng, rp, inv, y, y2;
    mpz_inits(g, Ng, rp, inv, y, y2, NULL);
    int rc = MR_FOUND;
    mpz_gcd(g, r, N);
    mpz_divexact(Ng, N, g);
    mpz_powm(y, u, Ng, mod);
    if (mpz_cmp_ui(y, 1) != 0) { rc = MR_NONE; goto done; }
    mpz_divexact(rp, r, g);
    if (mpz_cmp_ui(Ng, 1) == 0) {
        mpz_set_ui(y, 1);
    } else {
        mpz_invert(inv, rp, Ng);                     /* gcd(r/g, N/g) = 1 */
        mpz_powm(y, u, inv, mod);
    }
    if (mpz_cmp_ui(g, 1) != 0) {
        /* g | N, so its prime factors are the q the AMM step takes. */
        Expr* fargs[1] = { expr_bigint_normalize(expr_new_bigint_from_mpz(g)) };
        Expr* fact = internal_factorinteger(fargs, 1);
        if (!fact || !head_is(fact, SYM_List)) {
            if (fact) expr_free(fact);
            rc = MR_UNSUPPORTED; goto done;
        }
        for (size_t i = 0; i < fact->data.function.arg_count && rc == MR_FOUND; i++) {
            Expr* pr = fact->data.function.args[i];
            if (!head_is(pr, SYM_List) || pr->data.function.arg_count != 2 ||
                pr->data.function.args[0]->type != EXPR_INTEGER ||
                pr->data.function.args[1]->type != EXPR_INTEGER ||
                pr->data.function.args[0]->data.integer < 2) {
                rc = MR_UNSUPPORTED; break;
            }
            unsigned long q = (unsigned long)pr->data.function.args[0]->data.integer;
            for (int64_t k = 0; k < pr->data.function.args[1]->data.integer && rc == MR_FOUND; k++) {
                rc = amm_root(y2, y, q, N, mod, z_hint);
                mpz_set(y, y2);
            }
        }
        expr_free(fact);
    }
    if (rc == MR_FOUND) mpz_set(x, y);
done:
    mpz_clears(g, Ng, rp, inv, y, y2, NULL);
    return rc;
}

/* y with y^r ≡ u (mod p^k), u a unit mod p, k >= 1.  See step 3 above. */
static int unit_root_mod_pk(mpz_t y, const mpz_t u, const mpz_t r,
                            const mpz_t p, unsigned long k) {
    mpz_t pk, um, N, inv;
    mpz_inits(pk, um, N, inv, NULL);
    mpz_pow_ui(pk, p, k);
    mpz_mod(um, u, pk);
    int rc = MR_FOUND;
    if (mpz_cmp_ui(p, 2) != 0 && mpz_cmp_ui(r, 2) == 0) {
        /* Square root, odd p: Tonelli-Shanks mod p, Hensel-lifted (valid,
         * p ∤ 2 x0).  Kept for the established choice of root. */
        mpz_t up, x0;
        mpz_inits(up, x0, NULL);
        mpz_mod(up, um, p);
        if (!tonelli_shanks(x0, up, p)) rc = MR_NONE;
        else if (!hensel_lift(y, x0, r, um, p, k)) rc = MR_UNSUPPORTED;
        mpz_clears(up, x0, NULL);
    } else if (mpz_cmp_ui(p, 2) != 0) {
        mpz_pow_ui(N, p, k - 1);                     /* cyclic, order p^(k-1)(p-1) */
        mpz_sub_ui(inv, p, 1);
        mpz_mul(N, N, inv);
        rc = cyc_root(y, um, r, N, pk, NULL);
    } else if (k == 1) {
        mpz_set_ui(y, 1);
    } else if (mpz_odd_p(r)) {
        if (k == 2) {
            mpz_set(y, um);                          /* x^r = x on {1, 3} */
        } else {
            mpz_ui_pow_ui(N, 2, k - 2);              /* the group's exponent */
            mpz_invert(inv, r, N);
            mpz_powm(y, um, inv, pk);
        }
    } else if (mpz_fdiv_ui(um, 4) != 1) {
        rc = MR_NONE;                                /* an even power is 1 mod 4 */
    } else if (k == 2) {
        mpz_set_ui(y, 1);
    } else {
        mpz_t five;
        mpz_init_set_ui(five, 5);                    /* generates <5>, order 2^(k-2) */
        mpz_ui_pow_ui(N, 2, k - 2);
        rc = cyc_root(y, um, r, N, pk, five);
        mpz_clear(five);
    }
    mpz_clears(pk, um, N, inv, NULL);
    return rc;
}

/* x with x^r ≡ c (mod p^e), p prime; pe_out = p^e.  See step 2 above. */
static int rth_root_mod_pe(mpz_t x, mpz_t pe_out, const mpz_t c,
                           const mpz_t r, const mpz_t p, unsigned long e) {
    mpz_pow_ui(pe_out, p, e);
    /* A small prime power is settled exhaustively (its smallest root, or
     * none), as it always was. */
    if (mpz_cmp_ui(pe_out, MODROOT_BRUTE_LIMIT) <= 0)
        return modroot_brute(x, c, r, pe_out) ? MR_FOUND : MR_NONE;
    mpz_t cm, u, y, pw;
    mpz_inits(cm, u, y, pw, NULL);
    mpz_mod(cm, c, pe_out);
    int rc = MR_FOUND;
    if (mpz_sgn(cm) == 0) {
        mpz_set_ui(x, 0);
        goto done;
    }
    unsigned long v = (unsigned long)mpz_remove(u, cm, p);   /* v < e */
    unsigned long w = 0;
    if (v > 0) {
        if (!mpz_fits_ulong_p(r) || v % mpz_get_ui(r) != 0) { rc = MR_NONE; goto done; }
        w = v / mpz_get_ui(r);
    }
    rc = unit_root_mod_pk(y, u, r, p, e - v);
    if (rc != MR_FOUND) goto done;
    mpz_pow_ui(pw, p, w);
    mpz_mul(x, pw, y);
    mpz_mod(x, x, pe_out);
done:
    mpz_clears(cm, u, y, pw, NULL);
    return rc;
}

/* Top-level: x in [0, m) with x^r ≡ c (mod m), as MR_FOUND / MR_NONE (no
 * solution, proven) / MR_UNSUPPORTED (undecided: leave PowerMod alone). */
static int modular_root(mpz_t x, const mpz_t c, const mpz_t r, const mpz_t m) {
    if (mpz_cmp_ui(m, 1) == 0) { mpz_set_ui(x, 0); return MR_FOUND; }
    if (mpz_cmp_ui(m, 1) < 0)  return MR_UNSUPPORTED;
    if (mpz_sgn(r) <= 0) return MR_UNSUPPORTED;
    if (mpz_cmp_ui(r, 1) == 0) { mpz_mod(x, c, m); return MR_FOUND; }

    Expr* fact = NULL;
    if (!mpz_probab_prime_p(m, 25)) {
        /* internal_factorinteger takes ownership of the args array contents. */
        Expr* args[1] = { expr_bigint_normalize(expr_new_bigint_from_mpz(m)) };
        fact = internal_factorinteger(args, 1);
        if (!fact || !head_is(fact, SYM_List)) {
            if (fact) expr_free(fact);
            return MR_UNSUPPORTED;
        }
    }

    size_t n = fact ? fact->data.function.arg_count : 1;
    mpz_t x_acc, m_acc, pi, pe, xi, inv, diff;
    mpz_inits(x_acc, m_acc, pi, pe, xi, inv, diff, NULL);
    mpz_set_ui(x_acc, 0);
    mpz_set_ui(m_acc, 1);

    int rc = MR_FOUND;
    for (size_t i = 0; i < n && rc != MR_NONE; i++) {
        unsigned long ei = 1;
        if (fact) {
            Expr* pair = fact->data.function.args[i];
            if (!head_is(pair, SYM_List) || pair->data.function.arg_count != 2 ||
                !expr_is_integer_like(pair->data.function.args[0]) ||
                pair->data.function.args[1]->type != EXPR_INTEGER ||
                pair->data.function.args[1]->data.integer <= 0) {
                rc = MR_UNSUPPORTED; break;
            }
            mpz_t tp;
            expr_to_mpz(pair->data.function.args[0], tp);
            mpz_set(pi, tp);
            mpz_clear(tp);
            ei = (unsigned long)pair->data.function.args[1]->data.integer;
        } else {
            mpz_set(pi, m);
        }
        int ri = rth_root_mod_pe(xi, pe, c, r, pi, ei);
        if (ri == MR_NONE) { rc = MR_NONE; break; }
        if (ri != MR_FOUND) { rc = MR_UNSUPPORTED; continue; }  /* a later NONE still decides */
        if (rc != MR_FOUND) continue;
        /* CRT: x_acc (mod m_acc) and xi (mod pe), coprime moduli. */
        mpz_invert(inv, m_acc, pe);
        mpz_sub(diff, xi, x_acc);
        mpz_mul(diff, diff, inv);
        mpz_mod(diff, diff, pe);
        mpz_mul(diff, diff, m_acc);
        mpz_add(x_acc, x_acc, diff);
        mpz_mul(m_acc, m_acc, pe);
        mpz_mod(x_acc, x_acc, m_acc);
    }
    if (fact) expr_free(fact);

    if (rc != MR_NONE && mpz_fits_ulong_p(m) && mpz_get_ui(m) <= MODROOT_BRUTE_LIMIT) {
        /* A small modulus answers with its smallest root (and settles an
         * unsupported case either way). */
        if (modroot_brute(xi, c, r, m)) { mpz_set(x_acc, xi); rc = MR_FOUND; }
        else rc = MR_NONE;
    }
    if (rc == MR_FOUND) {
        /* Certify the root before returning it. */
        mpz_powm(xi, x_acc, r, m);
        mpz_mod(diff, c, m);
        if (mpz_cmp(xi, diff) != 0) rc = MR_UNSUPPORTED;
        else mpz_set(x, x_acc);
    }
    mpz_clears(x_acc, m_acc, pi, pe, xi, inv, diff, NULL);
    return rc;
}

/* PowerMod[a, b, m].  b may be an integer (incl. negative for the modular
 * inverse path) or Rational[p,q] (modular q-th root of a^p mod m). */
/* A number of any kind (what separates ::pmod from a symbolic call). */
static bool pm_is_number(Expr* e) {
    Expr *re, *im;
    return is_rational_like(e) || e->type == EXPR_REAL || is_complex(e, &re, &im)
#ifdef USE_MPFR
        || e->type == EXPR_MPFR
#endif
        ;
}

static bool pm_is_gaussian_int(Expr* e) {
    return expr_is_integer_like(e) || is_gaussian_integer(e);
}

static Expr* pm_emit_ninv(Expr* a, Expr* m) {
    char* sa = expr_to_string(a);
    char* sm = expr_to_string(m);
    mth_message_gated(g_arith_warnings_muted, "PowerMod", "ninv",
                      "%s is not invertible modulo %s.", sa ? sa : "?", sm ? sm : "?");
    free(sa); free(sm);
    return NULL;
}

/* PowerMod[a, b, m] over Z[i]: a Gaussian base or modulus, integer b. */
static Expr* pm_gaussian(Expr* a_expr, Expr* b_expr, Expr* m_expr) {
    Expr* fast;
    if (b_expr->type == EXPR_INTEGER && gi64_powmod(a_expr, b_expr->data.integer, m_expr, &fast))
        return fast;
    mpz_t ar, ai, mr, mi, b, rr, ri;
    if (!df_to_gaussian(a_expr, ar, ai)) return NULL;
    if (!df_to_gaussian(m_expr, mr, mi)) { mpz_clears(ar, ai, NULL); return NULL; }
    expr_to_mpz(b_expr, b);
    mpz_inits(rr, ri, NULL);
    Expr* out = NULL;
    if (gi_powmod(rr, ri, ar, ai, b, mr, mi)) {
        mpz_t one;
        mpz_init_set_ui(one, 1);
        out = gq_to_expr(rr, ri, one);
        mpz_clear(one);
    } else {
        pm_emit_ninv(a_expr, m_expr);
    }
    mpz_clears(ar, ai, mr, mi, b, rr, ri, NULL);
    return out;
}

/* PowerMod[a, b, m]: a^b mod m.
 *   integer b       square-and-multiply in GMP (mpz_powm); b < 0 inverts a
 *                   first, and a non-invertible a draws PowerMod::ninv;
 *   rational p/q    a q-th root of a^p mod m (PowerMod::root when none);
 *   Gaussian a or m the same over Z[i] for integer b, reducing with the
 *                   nearest-quotient Mod after every product.
 * As in the reference: an inexact argument, a non-integer base or modulus, or a
 * non-rational exponent draws PowerMod::pmod; a zero modulus PowerMod::divz; a
 * negative integer modulus leaves the call unevaluated; a wrong argument count
 * draws PowerMod::argr / ::argrx. Symbolic arguments stay unevaluated silently.
 * Listable, so the evaluator threads lists; an int64 buffer of bases takes the
 * ndint_powermod fast path. */
Expr* builtin_powermod(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 3) {
        if (argc == 1)
            mth_message_gated(g_arith_warnings_muted, "PowerMod", "argr",
                "PowerMod called with 1 argument; 3 arguments are expected.");
        else
            mth_message_gated(g_arith_warnings_muted, "PowerMod", "argrx",
                "PowerMod called with %zu arguments; 3 arguments are expected.", argc);
        return NULL;
    }
    Expr* a_expr = res->data.function.args[0];
    Expr* b_expr = res->data.function.args[1];
    Expr* m_expr = res->data.function.args[2];

    if (!is_ndarray(a_expr)) {
        /* Symbolic anywhere: silent. Numbers of the wrong kind: ::pmod. */
        if (!pm_is_number(a_expr) || !pm_is_number(b_expr) || !pm_is_number(m_expr))
            return NULL;
        if (!pm_is_gaussian_int(a_expr) || !pm_is_gaussian_int(m_expr) ||
            !is_rational_like(b_expr)) {
            char* call = expr_to_string(res);
            mth_message_gated(g_arith_warnings_muted, "PowerMod", "pmod",
                "The base and modulus arguments in %s must be Gaussian integers "
                "and the power must be an integer or rational.", call ? call : "?");
            free(call);
            return NULL;
        }
        bool m_zero = (m_expr->type == EXPR_INTEGER && m_expr->data.integer == 0);
        if (m_zero) {
            char* call = expr_to_string(res);
            mth_message_gated(g_arith_warnings_muted, "PowerMod", "divz",
                "The argument 0 in %s should be nonzero.", call ? call : "?");
            free(call);
            return NULL;
        }
        /* A negative integer modulus stays unevaluated. */
        if ((m_expr->type == EXPR_INTEGER && m_expr->data.integer < 0) ||
            (m_expr->type == EXPR_BIGINT && mpz_sgn(m_expr->data.bigint) < 0))
            return NULL;
        if (is_gaussian_integer(a_expr) || is_gaussian_integer(m_expr))
            return expr_is_integer_like(b_expr) ? pm_gaussian(a_expr, b_expr, m_expr) : NULL;
    }

    /* Buffer of bases, scalar exponent and modulus: square-and-multiply in
     * int64, no mpz_init per element. Anything it declines -- a negative
     * exponent, a modulus past the no-overflow bound, a non-integer buffer --
     * degrades explicitly, since expr_is_integer_like below is false for an
     * NDArray and would leave the call unevaluated. */
    if (is_ndarray(a_expr)) {
        Expr* nd = ndint_powermod(res);
        return nd ? nd : ndarray_delist_and_reeval(res);
    }

    if (!expr_is_integer_like(a_expr) || !expr_is_integer_like(m_expr)) return NULL;

    /* Machine base, non-negative exponent, modulus below 2^32: every product
     * of two residues fits uint64, so no GMP. */
    if (a_expr->type == EXPR_INTEGER && b_expr->type == EXPR_INTEGER &&
        m_expr->type == EXPR_INTEGER && b_expr->data.integer >= 0 &&
        m_expr->data.integer > 0 && m_expr->data.integer < ((int64_t)1 << 32)) {
        uint64_t mm = (uint64_t)m_expr->data.integer;
        int64_t av = a_expr->data.integer % (int64_t)mm;
        uint64_t base = (uint64_t)(av < 0 ? av + (int64_t)mm : av);
        uint64_t r = 1 % mm, e = (uint64_t)b_expr->data.integer;
        while (e) {
            if (e & 1) r = r * base % mm;
            e >>= 1;
            if (e) base = base * base % mm;
        }
        return expr_new_integer((int64_t)r);
    }

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
                return pm_emit_ninv(a_expr, m_expr);
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
        b_expr->data.function.head->data.symbol.name == SYM_Rational &&
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
                return pm_emit_ninv(a_expr, m_expr);
            }
            mpz_neg(p_neg, p);
            mpz_powm(c, inv, p_neg, m);
            mpz_clears(inv, p_neg, NULL);
        } else {
            mpz_powm(c, a, p, m);
        }

        int mr = modular_root(root, c, q, m);
        if (mr == MR_UNSUPPORTED) {
            /* Undecided (a root-degree prime past the discrete-log bound, or
             * an unfactorable modulus): no claim either way. */
            mpz_clears(p, q, c, root, a, m, NULL);
            return NULL;
        }
        if (mr == MR_NONE) {
            /* PowerMod::root -- x^q == a^p (mod m) has no solution. */
            char* sq = mpz_get_str(NULL, 10, q);
            char* sc = mpz_get_str(NULL, 10, c);
            char* sm = mpz_get_str(NULL, 10, m);
            mth_message_gated(g_arith_warnings_muted, "PowerMod", "root",
                "The equation x^%s = %s (mod %s) has no integer solutions.", sq, sc, sm);
            void (*gmp_free)(void*, size_t);
            mp_get_memory_functions(NULL, NULL, &gmp_free);
            gmp_free(sq, strlen(sq) + 1); gmp_free(sc, strlen(sc) + 1); gmp_free(sm, strlen(sm) + 1);
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

/* nt_gaussian.c -- shared Gaussian-integer & divisor enumeration (is_gaussian_integer, df_*, divisors_ordinary, divisors_gaussian).
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

/* True when e is an exact Gaussian integer Complex[a, b] with both real and
 * imaginary parts integer-like (machine int or BigInt). */
bool is_gaussian_integer(Expr* e) {
    Expr *re, *im;
    if (!is_complex(e, &re, &im)) return false;
    return expr_is_integer_like(re) && expr_is_integer_like(im);
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
bool df_factor_mpz(const mpz_t m, mpz_t** out_primes,
                          unsigned long** out_exps, size_t* out_n) {
    /* internal_factorinteger takes ownership of the args array contents, so the
     * constructed integer must NOT be freed here. */
    Expr* args[1] = { expr_bigint_normalize(expr_new_bigint_from_mpz(m)) };
    Expr* fact = internal_factorinteger(args, 1);
    if (!fact || fact->type != EXPR_FUNCTION ||
        fact->data.function.head->type != EXPR_SYMBOL ||
        fact->data.function.head->data.symbol.name != SYM_List) {
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
Expr* divisors_ordinary(const mpz_t nabs) {
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
void df_normalize_quadrant(mpz_t re, mpz_t im) {
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
bool df_to_gaussian(Expr* e, mpz_t a, mpz_t b) {
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

/* Factor the Gaussian integer z = a_in + b_in i into Gaussian prime powers
 * (gu[i] + gv[i] i)^ge[i].  The primes are NOT associate-normalised (a split
 * prime's conjugate appears as c - d i).  On success returns true and sets
 * *out_gu, *out_gv (arrays of init'd mpz_t), *out_ge and *out_gn; the caller
 * clears each prime and frees all three arrays.  A unit (N(z) == 1) factors as
 * the empty product (*out_gn == 0).  Returns false — allocating nothing — for
 * z == 0 or any factoring failure. */
bool df_gaussian_prime_factor(const mpz_t a_in, const mpz_t b_in,
                                     mpz_t** out_gu, mpz_t** out_gv,
                                     unsigned long** out_ge, size_t* out_gn) {
    mpz_t N, t;
    mpz_inits(N, t, NULL);
    mpz_mul(N, a_in, a_in);
    mpz_mul(t, b_in, b_in);
    mpz_add(N, N, t);
    mpz_clear(t);

    if (mpz_sgn(N) == 0) { mpz_clear(N); return false; }   /* z = 0 */

    mpz_t* primes;
    unsigned long* exps;
    size_t np;
    bool fok = df_factor_mpz(N, &primes, &exps, &np);
    mpz_clear(N);
    if (!fok) return false;

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
        return false;
    }
    *out_gu = gu;
    *out_gv = gv;
    *out_ge = ge;
    *out_gn = gn;
    return true;
}

/* Build the list of Gaussian divisors of z = a_in + b_in i, one representative
 * per associate class (first quadrant), sorted by (Re, Im).  Factors z into
 * Gaussian primes, enumerates the divisor lattice, and normalises each product
 * to its first-quadrant associate.  Returns a List Expr (just {1} for a unit),
 * or NULL on z == 0 / factoring failure. */
Expr* divisors_gaussian(const mpz_t a_in, const mpz_t b_in) {
    mpz_t* gu;
    mpz_t* gv;
    unsigned long* ge;
    size_t gn;
    if (!df_gaussian_prime_factor(a_in, b_in, &gu, &gv, &ge, &gn)) return NULL;

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

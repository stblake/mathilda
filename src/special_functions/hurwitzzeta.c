/* Mathilda -- the Hurwitz zeta function.
 *
 *   HurwitzZeta[s, a]   zeta(s,a) = Sum_{k>=0} (k+a)^-s        (Re s > 1)
 *
 * defined elsewhere by analytic continuation. HurwitzZeta agrees with the
 * two-argument Zeta for Re(a) > 0, but unlike Zeta it sums the *principal
 * branch* powers (k+a)^-s rather than ((k+a)^2)^(-s/2). The consequences:
 *
 *   - the two functions disagree for non-positive real a, and
 *   - HurwitzZeta retains the singular summands that Zeta discards, so it has
 *     poles at a = 0, -1, -2, ... .
 *
 * The evaluator routes each kind of argument to the cheapest exact or fastest
 * numeric path:
 *
 *   s == 1 (exact)             ->  ComplexInfinity (pole, for any a)
 *   a == 1                     ->  Zeta[s]          (Riemann closed forms)
 *   a == 1/2                   ->  (2^s - 1) Zeta[s]
 *   a positive integer m >= 2  ->  Zeta[s] - Sum_{k=1}^{m-1} k^-s
 *   a non-positive integer:
 *       s positive integer     ->  ComplexInfinity (pole)
 *       s non-positive integer ->  -BernoulliB[1-s, a]/(1-s)   (polynomial)
 *   any inexact operand        ->  Euler-Maclaurin complex-MPFR kernel
 *   everything else            ->  stays symbolic (return NULL)
 *
 * MPFR has no Hurwitz zeta, so the numeric kernel is implemented here from the
 * Euler-Maclaurin summation formula (DLMF 25.11.5):
 *
 *   zeta(s,a) = Sum_{k=0}^{N-1} (a+k)^-s
 *             + (a+N)^(1-s)/(s-1)
 *             + 1/2 (a+N)^-s
 *             + Sum_{j>=1} B_{2j}/(2j)! (s)_{2j-1} (a+N)^(-s-2j+1)
 *
 * with (s)_{2j-1} the rising factorial. N grows with the working precision and
 * |s|; the correction series is truncated at its optimal (smallest) term. The
 * kernel uses the principal branch for every (a+k)^-s, which is exactly the
 * HurwitzZeta convention. (The structure mirrors src/special_functions/zeta.c;
 * the self-contained Bernoulli cache and complex-MPFR toolkit are replicated
 * here so the two files stay independent.)
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "hurwitzzeta.h"
#include "sym_names.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "arithmetic.h"   /* is_rational, make_rational, is_complex, make_complex */
#include "numeric.h"      /* numeric_min_inexact_bits */
#include "attr.h"
#include "eval.h"          /* eval_and_free */
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#define HZRND MPFR_RNDN
#endif

/* M_PI is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Largest positive integer a for which exact Hurwitz expands to a partial sum. */
#define HZ_HURWITZ_A_CAP 100000L

/* ------------------------------------------------------------------ */
/* Small predicates / coercions                                        */
/* ------------------------------------------------------------------ */

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool hz_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if any leaf of `e` (including the parts of a Complex) is inexact. */
static bool hz_has_inexact(const Expr* e) {
    if (!e) return false;
    if (hz_is_inexact(e)) return true;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) return hz_is_inexact(re) || hz_is_inexact(im);
    return false;
}

/* Extract an exact machine integer from `e` (Integer, or BigInt fitting a
 * signed long). Inexact reals do NOT count. */
static bool hz_exact_int(const Expr* e, long* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (long)e->data.integer; return true; }
    if (e->type == EXPR_BIGINT && mpz_fits_slong_p(e->data.bigint)) {
        *out = mpz_get_si(e->data.bigint); return true;
    }
    return false;
}

/* True if `e` is exactly the rational 1/2. */
static bool hz_is_half(const Expr* e) {
    int64_t n, d;
    return is_rational(e, &n, &d) && n == 1 && d == 2;
}

/* ------------------------------------------------------------------ */
/* Small Expr builders (Zeta[s], etc.)                                 */
/* ------------------------------------------------------------------ */

static Expr* hz_zeta_of(const Expr* s) {
    return expr_new_function(expr_new_symbol(SYM_Zeta),
                             (Expr*[]){ expr_copy((Expr*)s) }, 1);
}

/* ------------------------------------------------------------------ */
/* Exact paths                                                         */
/* ------------------------------------------------------------------ */

/* HurwitzZeta[s, 1/2] = (2^s - 1) Zeta[s]. */
static Expr* hz_half(const Expr* s) {
    Expr* twos = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_new_integer(2), expr_copy((Expr*)s) }, 2);
    Expr* coeff = expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ twos, expr_new_integer(-1) }, 2);
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ coeff, hz_zeta_of(s) }, 2);
    return eval_and_free(prod);
}

/* Exact Hurwitz at a positive integer a = m:
 *   zeta(s, m) = Zeta[s] - Sum_{k=1}^{m-1} k^-s.
 * Valid for symbolic / exact / integer s. Returns NULL past the cap. */
static Expr* hz_int_a(const Expr* s, long m) {
    if (m < 1 || m > HZ_HURWITZ_A_CAP) return NULL;

    Expr* zs = hz_zeta_of(s);
    if (m == 1) return eval_and_free(zs);          /* zeta(s, 1) = zeta(s) */

    size_t terms = (size_t)m;                       /* zs + (m-1) subtracted terms */
    Expr** parts = (Expr**)malloc(terms * sizeof(Expr*));
    if (!parts) { expr_free(zs); return NULL; }
    parts[0] = zs;
    for (long k = 1; k < m; k++) {
        /* -(k^-s) = Times[-1, Power[k, Times[-1, s]]]. */
        Expr* negs = expr_new_function(expr_new_symbol(SYM_Times),
                         (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)s) }, 2);
        Expr* kpow = expr_new_function(expr_new_symbol(SYM_Power),
                         (Expr*[]){ expr_new_integer(k), negs }, 2);
        parts[(size_t)k] = expr_new_function(expr_new_symbol(SYM_Times),
                         (Expr*[]){ expr_new_integer(-1), kpow }, 2);
    }
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), parts, terms);
    free(parts);
    return eval_and_free(sum);
}

/* HurwitzZeta[-m, a] = -BernoulliB[m+1, a]/(m+1) for non-positive integer s
 * (s = -m, m >= 0). Entire in a, so finite even at non-positive integer a
 * (e.g. HurwitzZeta[0, 0] = 1/2). `sn` is s as a machine integer (<= 0). */
static Expr* hz_neg_int_s(long sn, const Expr* a) {
    long k = 1 - sn;                                /* m + 1 = 1 - s, >= 1 */
    Expr* bern = expr_new_function(expr_new_symbol(SYM_BernoulliB),
                     (Expr*[]){ expr_new_integer(k), expr_copy((Expr*)a) }, 2);
    Expr* coeff = make_rational(-1, k);             /* -1/(m+1) */
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ coeff, bern }, 2);
    return eval_and_free(prod);
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ */
/* Bernoulli numbers B_k (exact, lazily cached, process-lifetime)      */
/*                                                                     */
/* Recurrence: B_0 = 1, and for m >= 1                                 */
/*   B_m = -1/(m+1) Sum_{k=0}^{m-1} C(m+1, k) B_k.                      */
/* Used only by the Euler-Maclaurin correction series.                 */
/* ------------------------------------------------------------------ */

static mpq_t* g_hzbern = NULL;
static size_t g_hzbern_len = 0;

static void hzbern_ensure(size_t upto) {
    if (upto + 1 <= g_hzbern_len) return;
    size_t newlen = upto + 1;
    mpq_t* grown = (mpq_t*)realloc(g_hzbern, newlen * sizeof(mpq_t));
    if (!grown) return; /* out of memory: leave cache as-is, callers degrade */
    g_hzbern = grown;

    mpq_t sum, term, factor;
    mpq_inits(sum, term, factor, (mpq_ptr)0);
    mpz_t binz;
    mpz_init(binz);

    for (size_t m = g_hzbern_len; m < newlen; m++) {
        mpq_init(g_hzbern[m]);
        if (m == 0) { mpq_set_ui(g_hzbern[m], 1, 1); continue; }
        mpq_set_ui(sum, 0, 1);
        for (size_t k = 0; k < m; k++) {
            mpz_bin_uiui(binz, (unsigned long)(m + 1), (unsigned long)k);
            mpq_set_z(factor, binz);
            mpq_mul(term, factor, g_hzbern[k]);
            mpq_add(sum, sum, term);
        }
        mpq_set_si(term, -1, (unsigned long)(m + 1));
        mpq_canonicalize(term);
        mpq_mul(g_hzbern[m], sum, term);
        mpq_canonicalize(g_hzbern[m]);
    }
    g_hzbern_len = newlen;
    mpz_clear(binz);
    mpq_clears(sum, term, factor, (mpq_ptr)0);
}

/* B_idx into `out`. Odd indices above 1 are exactly zero. */
static void hzbern_get_q(mpq_t out, size_t idx) {
    if (idx > 1 && (idx & 1u)) { mpq_set_ui(out, 0, 1); return; }
    hzbern_ensure(idx);
    if (idx < g_hzbern_len) mpq_set(out, g_hzbern[idx]);
    else mpq_set_ui(out, 0, 1);
}

/* (2j)! into out. */
static void hz_factorial(mpz_t out, long n) {
    mpz_set_ui(out, 1);
    for (long i = 2; i <= n; i++) mpz_mul_ui(out, out, (unsigned long)i);
}

/* ------------------------------------------------------------------ */
/* Complex-MPFR toolkit (pairs of mpfr_t; no MPC library available).   */
/* Mirrors the `zcx` helpers in zeta.c -- alias-safe, explicit         */
/* working precision `p`.                                              */
/* ------------------------------------------------------------------ */

typedef struct { mpfr_t re, im; } hcx;

static void hcx_init(hcx* z, mpfr_prec_t p) { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void hcx_clear(hcx* z)               { mpfr_clear(z->re);    mpfr_clear(z->im);    }
static void hcx_set(hcx* d, const hcx* s)   { mpfr_set(d->re, s->re, HZRND); mpfr_set(d->im, s->im, HZRND); }

static void hcx_add(hcx* out, const hcx* a, const hcx* b) {
    mpfr_add(out->re, a->re, b->re, HZRND);
    mpfr_add(out->im, a->im, b->im, HZRND);
}

static void hcx_abs(mpfr_t mag, const hcx* z) { mpfr_hypot(mag, z->re, z->im, HZRND); }

static void hcx_mul(hcx* out, const hcx* a, const hcx* b, mpfr_prec_t p) {
    mpfr_t ac, bd, ad, bc;
    mpfr_inits2(p, ac, bd, ad, bc, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, HZRND);
    mpfr_mul(bd, a->im, b->im, HZRND);
    mpfr_mul(ad, a->re, b->im, HZRND);
    mpfr_mul(bc, a->im, b->re, HZRND);
    mpfr_sub(out->re, ac, bd, HZRND);
    mpfr_add(out->im, ad, bc, HZRND);
    mpfr_clears(ac, bd, ad, bc, (mpfr_ptr)0);
}

/* out = a * r, r a real scalar. */
static void hcx_mul_r(hcx* out, const hcx* a, const mpfr_t r) {
    mpfr_mul(out->re, a->re, r, HZRND);
    mpfr_mul(out->im, a->im, r, HZRND);
}

static void hcx_div(hcx* out, const hcx* a, const hcx* b, mpfr_prec_t p) {
    mpfr_t den, t1, t2, nr, ni;
    mpfr_inits2(p, den, t1, t2, nr, ni, (mpfr_ptr)0);
    mpfr_mul(t1, b->re, b->re, HZRND);
    mpfr_mul(t2, b->im, b->im, HZRND);
    mpfr_add(den, t1, t2, HZRND);                 /* |b|^2 */
    mpfr_mul(t1, a->re, b->re, HZRND);
    mpfr_mul(t2, a->im, b->im, HZRND);
    mpfr_add(nr, t1, t2, HZRND);                  /* ac + bd */
    mpfr_mul(t1, a->im, b->re, HZRND);
    mpfr_mul(t2, a->re, b->im, HZRND);
    mpfr_sub(ni, t1, t2, HZRND);                  /* bc - ad */
    mpfr_div(out->re, nr, den, HZRND);
    mpfr_div(out->im, ni, den, HZRND);
    mpfr_clears(den, t1, t2, nr, ni, (mpfr_ptr)0);
}

static void hcx_exp(hcx* out, const hcx* z, mpfr_prec_t p) {
    mpfr_t ea, c, s;
    mpfr_inits2(p, ea, c, s, (mpfr_ptr)0);
    mpfr_exp(ea, z->re, HZRND);
    mpfr_sin_cos(s, c, z->im, HZRND);
    mpfr_mul(out->re, ea, c, HZRND);
    mpfr_mul(out->im, ea, s, HZRND);
    mpfr_clears(ea, c, s, (mpfr_ptr)0);
}

/* out = Log(z), principal branch. */
static void hcx_log(hcx* out, const hcx* z, mpfr_prec_t p) {
    mpfr_t mag, ang;
    mpfr_inits2(p, mag, ang, (mpfr_ptr)0);
    mpfr_hypot(mag, z->re, z->im, HZRND);
    mpfr_atan2(ang, z->im, z->re, HZRND);
    mpfr_log(out->re, mag, HZRND);
    mpfr_set(out->im, ang, HZRND);
    mpfr_clears(mag, ang, (mpfr_ptr)0);
}

/* out = base^expo (principal branch) = exp(expo * Log(base)). */
static void hcx_pow(hcx* out, const hcx* base, const hcx* expo, mpfr_prec_t p) {
    hcx lg, prod;
    hcx_init(&lg, p); hcx_init(&prod, p);
    hcx_log(&lg, base, p);
    hcx_mul(&prod, expo, &lg, p);
    hcx_exp(out, &prod, p);
    hcx_clear(&lg); hcx_clear(&prod);
}

/* ------------------------------------------------------------------ */
/* Euler-Maclaurin Hurwitz zeta kernel                                 */
/* ------------------------------------------------------------------ */

/* Result codes from the numeric kernel. */
enum { HZ_OK = 0, HZ_NONFINITE = 1, HZ_SINGULAR = 2 };

/* Number of head terms N for working precision wp and |s|. */
static unsigned long hz_em_terms(mpfr_prec_t wp, const hcx* s) {
    double digits = (double)wp * 0.30103 + 5.0;
    double smag = hypot(mpfr_get_d(s->re, HZRND), mpfr_get_d(s->im, HZRND));
    double n = digits + smag + 10.0;
    if (n < 8.0) n = 8.0;
    if (n > 200000.0) n = 200000.0;
    return (unsigned long)n;
}

/* zeta(s, a) for complex s, a into `out`. Returns HZ_OK on success,
 * HZ_NONFINITE if the result is non-finite (e.g. the s = 1 pole), or
 * HZ_SINGULAR if a head term (k + a) hits exactly zero while Re(s) > 0 (a pole
 * at a non-positive integer a). */
static int hz_hurwitz_cx(hcx* out, const hcx* s, const hcx* a, mpfr_prec_t wp) {
    unsigned long N = hz_em_terms(wp, s);
    bool singular = false;

    hcx S, negs, base, term, w, P0, f, winv2, w2, sm1, r, tmp, tmp2;
    hcx_init(&S, wp); hcx_init(&negs, wp); hcx_init(&base, wp); hcx_init(&term, wp);
    hcx_init(&w, wp); hcx_init(&P0, wp); hcx_init(&f, wp); hcx_init(&winv2, wp);
    hcx_init(&w2, wp); hcx_init(&sm1, wp); hcx_init(&r, wp); hcx_init(&tmp, wp);
    hcx_init(&tmp2, wp);

    mpfr_t eps, mag, prevmag, sabs, cj, bq_d;
    mpfr_inits2(wp, eps, mag, prevmag, sabs, cj, bq_d, (mpfr_ptr)0);
    mpfr_set_ui(eps, 1, HZRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 12 ? wp - 8 : 1), HZRND);

    /* -s */
    mpfr_neg(negs.re, s->re, HZRND);
    mpfr_neg(negs.im, s->im, HZRND);

    /* Head sum S = Sum_{k=0}^{N-1} (a+k)^-s. A zero base is a pole when
     * Re(s) > 0; for Re(s) <= 0 the principal-branch power is 0, so skip it. */
    mpfr_set_ui(S.re, 0, HZRND); mpfr_set_ui(S.im, 0, HZRND);
    for (unsigned long k = 0; k < N; k++) {
        mpfr_add_ui(base.re, a->re, k, HZRND);
        mpfr_set(base.im, a->im, HZRND);
        hcx_abs(mag, &base);
        if (mpfr_sgn(mag) == 0) {                  /* k + a == 0 */
            if (mpfr_sgn(s->re) > 0) singular = true;
            continue;
        }
        hcx_pow(&term, &base, &negs, wp);
        hcx_add(&S, &S, &term);
    }

    /* w = a + N. */
    mpfr_add_ui(w.re, a->re, N, HZRND);
    mpfr_set(w.im, a->im, HZRND);

    /* P0 = w^-s. */
    hcx_pow(&P0, &w, &negs, wp);

    /* + w^(1-s)/(s-1) = (P0 * w)/(s-1). */
    mpfr_sub_ui(sm1.re, s->re, 1, HZRND);
    mpfr_set(sm1.im, s->im, HZRND);
    hcx_mul(&tmp, &P0, &w, wp);
    hcx_div(&tmp2, &tmp, &sm1, wp);
    hcx_add(&S, &S, &tmp2);

    /* + 1/2 w^-s. */
    mpfr_set_d(cj, 0.5, HZRND);
    hcx_mul_r(&tmp, &P0, cj);
    hcx_add(&S, &S, &tmp);

    /* Correction series.  f_j = w^-(2j-1), starting f_1 = 1/w; advance by
     * winv2 = 1/w^2 each step.  r_j = (s)_{2j-1}, starting r_1 = s. */
    {
        hcx one; hcx_init(&one, wp);
        mpfr_set_ui(one.re, 1, HZRND); mpfr_set_ui(one.im, 0, HZRND);
        hcx_div(&f, &one, &w, wp);                 /* 1/w */
        hcx_mul(&w2, &w, &w, wp);
        hcx_div(&winv2, &one, &w2, wp);            /* 1/w^2 */
        hcx_clear(&one);
    }
    hcx_set(&r, s);                                /* (s)_1 = s */

    mpfr_set_inf(prevmag, 1);
    mpq_t bq;
    mpq_init(bq);
    mpz_t fac;
    mpz_init(fac);

    unsigned long jmax = (unsigned long)wp + 8;
    for (unsigned long j = 1; j <= jmax; j++) {
        /* coeff c_j = B_{2j}/(2j)!. */
        hzbern_get_q(bq, (size_t)(2 * j));
        hz_factorial(fac, (long)(2 * j));
        mpfr_set_q(cj, bq, HZRND);
        mpfr_set_z(bq_d, fac, HZRND);
        mpfr_div(cj, cj, bq_d, HZRND);             /* B_{2j}/(2j)! */

        /* term = c_j * r * P0 * f. */
        hcx_mul(&tmp, &r, &P0, wp);
        hcx_mul(&tmp2, &tmp, &f, wp);
        hcx_mul_r(&term, &tmp2, cj);
        hcx_abs(mag, &term);

        /* Optimal truncation: stop once the terms start growing again. */
        if (j > 1 && mpfr_cmp(mag, prevmag) > 0) break;
        hcx_add(&S, &S, &term);
        mpfr_set(prevmag, mag, HZRND);

        /* Converged? |term| < eps * |S|. */
        hcx_abs(sabs, &S);
        mpfr_mul(sabs, sabs, eps, HZRND);
        if (mpfr_cmp(mag, sabs) < 0) break;

        /* r *= (s + 2j-1)(s + 2j). */
        mpfr_add_ui(tmp.re, s->re, 2 * j - 1, HZRND); mpfr_set(tmp.im, s->im, HZRND);
        hcx_mul(&tmp2, &r, &tmp, wp);
        mpfr_add_ui(tmp.re, s->re, 2 * j, HZRND);     mpfr_set(tmp.im, s->im, HZRND);
        hcx_mul(&r, &tmp2, &tmp, wp);
        /* f *= 1/w^2. */
        hcx_mul(&tmp, &f, &winv2, wp);
        hcx_set(&f, &tmp);
    }

    mpz_clear(fac);
    mpq_clear(bq);
    hcx_set(out, &S);
    bool finite = mpfr_number_p(out->re) && mpfr_number_p(out->im);

    mpfr_clears(eps, mag, prevmag, sabs, cj, bq_d, (mpfr_ptr)0);
    hcx_clear(&S); hcx_clear(&negs); hcx_clear(&base); hcx_clear(&term);
    hcx_clear(&w); hcx_clear(&P0); hcx_clear(&f); hcx_clear(&winv2);
    hcx_clear(&w2); hcx_clear(&sm1); hcx_clear(&r); hcx_clear(&tmp); hcx_clear(&tmp2);

    if (singular) return HZ_SINGULAR;
    return finite ? HZ_OK : HZ_NONFINITE;
}

/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool hz_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, HZRND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          HZRND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        HZRND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          HZRND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, HZRND);
        mpfr_div_si(out, out, (long)d, HZRND);
        return true;
    }
    return false;
}

/* Result precision (bits) under Mathematica contagion: the minimum precision
 * among the inexact leaves of s and a, floored at machine (53). */
static mpfr_prec_t hz_out_prec(const Expr* s, const Expr* a) {
    long ps = numeric_min_inexact_bits(s);
    long pa = numeric_min_inexact_bits(a);
    long m  = (ps && pa) ? (ps < pa ? ps : pa) : (ps ? ps : pa);
    return m < 53 ? 53 : (mpfr_prec_t)m;
}

/* Build the numeric result: Real / Complex[Real,Real] for machine precision,
 * MPFR parts for arbitrary precision.  `real_only` forces a real result. */
static Expr* hz_result(const mpfr_t re, const mpfr_t im,
                       mpfr_prec_t out_prec, bool real_only) {
    if (out_prec <= 53) {
        Expr* rr = expr_new_real(mpfr_get_d(re, HZRND));
        if (real_only) return rr;
        return make_complex(rr, expr_new_real(mpfr_get_d(im, HZRND)));
    }
    Expr* rr = expr_new_mpfr_bits(out_prec);
    mpfr_set(rr->data.mpfr, re, HZRND);
    if (real_only) return rr;
    Expr* ii = expr_new_mpfr_bits(out_prec);
    mpfr_set(ii->data.mpfr, im, HZRND);
    return make_complex(rr, ii);
}

/* Decompose a numeric value into borrowed (re, im). */
static void hz_decompose(Expr* v, Expr* zero, Expr** re, Expr** im) {
    Expr *r, *i;
    if (is_complex(v, &r, &i)) { *re = r; *im = i; }
    else { *re = v; *im = zero; }
}

/* Numeric HurwitzZeta[s, a]. At least one of s, a carries an inexact part.
 * Returns Real/Complex/MPFR, ComplexInfinity at a pole, or NULL if the parts
 * are not numeric. */
static Expr* hz_numeric(Expr* s, Expr* a) {
    Expr* zero = expr_new_integer(0);
    Expr *sre, *sim, *are, *aim;
    hz_decompose(s, zero, &sre, &sim);
    hz_decompose(a, zero, &are, &aim);

    mpfr_prec_t out_prec = hz_out_prec(s, a);
    mpfr_prec_t wp = out_prec + 64;

    bool real_only = (sim == zero) && (aim == zero);

    hcx sc, ac, res;
    hcx_init(&sc, wp); hcx_init(&ac, wp); hcx_init(&res, wp);
    bool ok = hz_set_mpfr(sc.re, sre) && hz_set_mpfr(sc.im, sim) &&
              hz_set_mpfr(ac.re, are) && hz_set_mpfr(ac.im, aim);
    Expr* out = NULL;
    if (ok) {
        int code = hz_hurwitz_cx(&res, &sc, &ac, wp);
        if (code == HZ_OK) out = hz_result(res.re, res.im, out_prec, real_only);
        else               out = expr_new_symbol(SYM_ComplexInfinity);
    }
    hcx_clear(&sc); hcx_clear(&ac); hcx_clear(&res);
    expr_free(zero);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* HurwitzZeta[s, a]                                                   */
/* ------------------------------------------------------------------ */

static Expr* hz_two_arg(Expr* s, Expr* a) {
    bool s_inexact = hz_has_inexact(s);
    bool a_inexact = hz_has_inexact(a);

    /* Exact reductions (only when neither operand is inexact). */
    if (!s_inexact && !a_inexact) {
        long sn, m;
        bool s_is_int = hz_exact_int(s, &sn);

        /* s == 1: pole for any a. */
        if (s_is_int && sn == 1) return expr_new_symbol(SYM_ComplexInfinity);

        /* a == 1/2: (2^s - 1) Zeta[s]. */
        if (hz_is_half(a)) return hz_half(s);

        if (hz_exact_int(a, &m)) {
            if (m >= 1) {                          /* positive integer a */
                Expr* closed = hz_int_a(s, m);
                if (closed) return closed;
            } else {                               /* non-positive integer a */
                if (s_is_int) {
                    if (sn >= 1) return expr_new_symbol(SYM_ComplexInfinity);
                    return hz_neg_int_s(sn, a);    /* finite Bernoulli value */
                }
                /* symbolic / non-integer s at a pole point: stay symbolic. */
            }
        }
        return NULL; /* leave symbolic */
    }

    /* a == 1 with inexact s: HurwitzZeta[s, 1] = Zeta[s]. */
    if (!a_inexact && a->type == EXPR_INTEGER && a->data.integer == 1)
        return eval_and_free(hz_zeta_of(s));

#ifdef USE_MPFR
    /* Numeric (real or complex) kernel: at least one inexact operand. */
    {
        Expr* out = hz_numeric(s, a);
        if (out) return out;
    }
#endif

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

static Expr* hz_emit_argrx(size_t argc) {
    fprintf(stderr,
            "HurwitzZeta::argrx: HurwitzZeta called with %zu argument%s; "
            "2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_hurwitzzeta(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 2) return hz_two_arg(args[0], args[1]);
    return hz_emit_argrx(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void hurwitzzeta_init(void) {
    symtab_add_builtin("HurwitzZeta", builtin_hurwitzzeta);
    symtab_get_def("HurwitzZeta")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

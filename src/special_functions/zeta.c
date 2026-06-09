/* Mathilda -- the Riemann and Hurwitz zeta functions.
 *
 *   Zeta[s]      Riemann zeta      zeta(s)   = Sum_{k>=1} k^-s          (Re s > 1)
 *   Zeta[s, a]   Hurwitz zeta      zeta(s,a) = Sum_{k>=0} (k+a)^-s      (Re s > 1)
 *
 * Both are defined elsewhere by analytic continuation; the evaluator routes
 * each kind of argument to the cheapest exact or fastest numeric path:
 *
 *   exact integer s         ->  closed form:
 *                                 s = 1        : ComplexInfinity (pole)
 *                                 s = 0        : -1/2
 *                                 s = 2n > 0   : rational * Pi^(2n)   (Bernoulli)
 *                                 s = -m < 0   : rational            (Bernoulli)
 *                                 s = 2n+1 > 0 : stays symbolic (no closed form)
 *   exact Hurwitz, integer a -> Zeta[s] - Sum_{k=1}^{a-1} k^-s
 *   machine / MPFR real s    -> MPFR mpfr_zeta (Riemann only)
 *   complex s, or any a != 1 -> Euler-Maclaurin complex-MPFR kernel
 *   everything else          -> stays symbolic (return NULL)
 *
 * MPFR provides mpfr_zeta for real Riemann zeta only -- it has no Hurwitz and
 * no complex zeta -- so the Hurwitz / complex kernel is implemented here from
 * the Euler-Maclaurin summation formula (DLMF 25.11.5):
 *
 *   zeta(s,a) = Sum_{k=0}^{N-1} (a+k)^-s
 *             + (a+N)^(1-s)/(s-1)
 *             + 1/2 (a+N)^-s
 *             + Sum_{j>=1} B_{2j}/(2j)! (s)_{2j-1} (a+N)^(-s-2j+1)
 *
 * with (s)_{2j-1} the rising factorial. N is chosen from the working precision
 * and |s|; the correction series is truncated at its optimal (smallest) term.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "zeta.h"

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
#define ZRND MPFR_RNDN
#endif

/* M_PI is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Largest |s| (exact integer) for which we build the Bernoulli closed form;
 * beyond this the Bernoulli numbers get unwieldy, so stay symbolic / numeric. */
#define ZETA_EXACT_INT_CAP 10000L
/* Largest positive integer a for which exact Hurwitz expands to a partial sum. */
#define ZETA_HURWITZ_A_CAP 100000L

/* ------------------------------------------------------------------ */
/* Small predicates / coercions                                        */
/* ------------------------------------------------------------------ */

/* True if `e` is exactly the symbol `name`. */
static bool zeta_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol, name) == 0;
}

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool zeta_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* Extract an exact machine integer from `e` (Integer, or BigInt fitting a
 * signed long). Inexact reals do NOT count -- the exact path is only for
 * genuinely exact integer arguments. */
static bool zeta_exact_int(const Expr* e, long* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (long)e->data.integer; return true; }
    if (e->type == EXPR_BIGINT && mpz_fits_slong_p(e->data.bigint)) {
        *out = mpz_get_si(e->data.bigint); return true;
    }
    return false;
}

/* ------------------------------------------------------------------ */
/* Bernoulli numbers B_k (exact, lazily cached, process-lifetime)      */
/*                                                                     */
/* Recurrence: B_0 = 1, and for m >= 1                                 */
/*   B_m = -1/(m+1) Sum_{k=0}^{m-1} C(m+1, k) B_k.                      */
/* (Same construction as src/polygamma.c; replicated to stay            */
/* self-contained -- the cache there is file-static.)                  */
/* ------------------------------------------------------------------ */

static mpq_t* g_zbern = NULL;
static size_t g_zbern_len = 0;

static void zbern_ensure(size_t upto) {
    if (upto + 1 <= g_zbern_len) return;
    size_t newlen = upto + 1;
    mpq_t* grown = (mpq_t*)realloc(g_zbern, newlen * sizeof(mpq_t));
    if (!grown) return; /* out of memory: leave cache as-is, callers degrade */
    g_zbern = grown;

    mpq_t sum, term, factor;
    mpq_inits(sum, term, factor, (mpq_ptr)0);
    mpz_t binz;
    mpz_init(binz);

    for (size_t m = g_zbern_len; m < newlen; m++) {
        mpq_init(g_zbern[m]);
        if (m == 0) { mpq_set_ui(g_zbern[m], 1, 1); continue; }
        mpq_set_ui(sum, 0, 1);
        for (size_t k = 0; k < m; k++) {
            mpz_bin_uiui(binz, (unsigned long)(m + 1), (unsigned long)k);
            mpq_set_z(factor, binz);
            mpq_mul(term, factor, g_zbern[k]);
            mpq_add(sum, sum, term);
        }
        mpq_set_si(term, -1, (unsigned long)(m + 1));
        mpq_canonicalize(term);
        mpq_mul(g_zbern[m], sum, term);
        mpq_canonicalize(g_zbern[m]);
    }
    g_zbern_len = newlen;
    mpz_clear(binz);
    mpq_clears(sum, term, factor, (mpq_ptr)0);
}

/* B_idx into `out`. Odd indices above 1 are exactly zero. */
static void zbern_get_q(mpq_t out, size_t idx) {
    if (idx > 1 && (idx & 1u)) { mpq_set_ui(out, 0, 1); return; }
    zbern_ensure(idx);
    if (idx < g_zbern_len) mpq_set(out, g_zbern[idx]);
    else mpq_set_ui(out, 0, 1);
}

/* ------------------------------------------------------------------ */
/* Exact-value builders                                                */
/* ------------------------------------------------------------------ */

/* (2j)! into out. */
static void zeta_factorial(mpz_t out, long n) {
    mpz_set_ui(out, 1);
    for (long i = 2; i <= n; i++) mpz_mul_ui(out, out, (unsigned long)i);
}

/* Build a canonical Integer/BigInt/Rational from an mpq via Times[num,
 * Power[den,-1]] so the evaluator normalises arbitrary-size components. */
static Expr* zeta_expr_from_mpq(const mpq_t q) {
    mpz_t num, den;
    mpz_init(num); mpz_init(den);
    mpq_get_num(num, q);
    mpq_get_den(den, q);
    Expr* en = expr_bigint_normalize(expr_new_bigint_from_mpz(num));
    Expr* out;
    if (mpz_cmp_ui(den, 1) == 0) {
        out = en;
    } else {
        Expr* ed = expr_bigint_normalize(expr_new_bigint_from_mpz(den));
        Expr* inv = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ ed, expr_new_integer(-1) }, 2);
        out = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ en, inv }, 2));
    }
    mpz_clear(num); mpz_clear(den);
    return out;
}

/* zeta(2n) = (-1)^(n+1) B_{2n} 2^(2n-1) Pi^(2n) / (2n)!,  for n >= 1.
 * Returns rational * Pi^(2n), or NULL if 2n exceeds the cap. */
static Expr* zeta_exact_even(long two_n) {
    long n = two_n / 2;
    if (two_n > ZETA_EXACT_INT_CAP) return NULL;

    mpq_t c, tmp;
    mpz_t pw, f2n;
    mpq_inits(c, tmp, (mpq_ptr)0);
    mpz_inits(pw, f2n, (mpz_ptr)0);

    zbern_get_q(c, (size_t)two_n);                 /* B_{2n} */
    mpz_ui_pow_ui(pw, 2, (unsigned long)(two_n - 1));   /* 2^(2n-1) */
    mpq_set_z(tmp, pw);
    mpq_mul(c, c, tmp);
    zeta_factorial(f2n, two_n);                    /* (2n)! */
    mpq_set_z(tmp, f2n);
    mpq_div(c, c, tmp);
    if ((n % 2) == 0) mpq_neg(c, c);               /* (-1)^(n+1) */
    mpq_canonicalize(c);

    Expr* crat = zeta_expr_from_mpq(c);
    Expr* pipow = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_symbol("Pi"), expr_new_integer(two_n) }, 2);
    Expr* out = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ crat, pipow }, 2));

    mpq_clears(c, tmp, (mpq_ptr)0);
    mpz_clears(pw, f2n, (mpz_ptr)0);
    return out;
}

/* zeta(-m) = (-1)^m B_{m+1}/(m+1),  for m >= 1 (rational; 0 at even m). */
static Expr* zeta_exact_negint(long m) {
    if (m > ZETA_EXACT_INT_CAP) return NULL;
    if ((m & 1L) == 0) return expr_new_integer(0);     /* even m: B_{odd}=0 */

    mpq_t b, out;
    mpq_inits(b, out, (mpq_ptr)0);
    zbern_get_q(b, (size_t)(m + 1));                   /* B_{m+1} */
    mpq_set_si(out, 1, (unsigned long)(m + 1));         /* 1/(m+1) */
    mpq_mul(out, out, b);
    if (m & 1L) mpq_neg(out, out);                      /* (-1)^m, m odd */
    mpq_canonicalize(out);
    Expr* r = zeta_expr_from_mpq(out);
    mpq_clears(b, out, (mpq_ptr)0);
    return r;
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ */
/* Complex-MPFR toolkit (pairs of mpfr_t; no MPC library available).   */
/* Mirrors the `gcx` helpers in src/gamma.c -- alias-safe, explicit    */
/* working precision `p`.                                              */
/* ------------------------------------------------------------------ */

typedef struct { mpfr_t re, im; } zcx;

static void zcx_init(zcx* z, mpfr_prec_t p) { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void zcx_clear(zcx* z)               { mpfr_clear(z->re);    mpfr_clear(z->im);    }
static void zcx_set(zcx* d, const zcx* s)   { mpfr_set(d->re, s->re, ZRND); mpfr_set(d->im, s->im, ZRND); }

static void zcx_add(zcx* out, const zcx* a, const zcx* b) {
    mpfr_add(out->re, a->re, b->re, ZRND);
    mpfr_add(out->im, a->im, b->im, ZRND);
}

static void zcx_abs(mpfr_t mag, const zcx* z) { mpfr_hypot(mag, z->re, z->im, ZRND); }

static void zcx_mul(zcx* out, const zcx* a, const zcx* b, mpfr_prec_t p) {
    mpfr_t ac, bd, ad, bc;
    mpfr_inits2(p, ac, bd, ad, bc, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, ZRND);
    mpfr_mul(bd, a->im, b->im, ZRND);
    mpfr_mul(ad, a->re, b->im, ZRND);
    mpfr_mul(bc, a->im, b->re, ZRND);
    mpfr_sub(out->re, ac, bd, ZRND);
    mpfr_add(out->im, ad, bc, ZRND);
    mpfr_clears(ac, bd, ad, bc, (mpfr_ptr)0);
}

/* out = a * r, r a real scalar. */
static void zcx_mul_r(zcx* out, const zcx* a, const mpfr_t r) {
    mpfr_mul(out->re, a->re, r, ZRND);
    mpfr_mul(out->im, a->im, r, ZRND);
}

static void zcx_div(zcx* out, const zcx* a, const zcx* b, mpfr_prec_t p) {
    mpfr_t den, t1, t2, nr, ni;
    mpfr_inits2(p, den, t1, t2, nr, ni, (mpfr_ptr)0);
    mpfr_mul(t1, b->re, b->re, ZRND);
    mpfr_mul(t2, b->im, b->im, ZRND);
    mpfr_add(den, t1, t2, ZRND);                 /* |b|^2 */
    mpfr_mul(t1, a->re, b->re, ZRND);
    mpfr_mul(t2, a->im, b->im, ZRND);
    mpfr_add(nr, t1, t2, ZRND);                  /* ac + bd */
    mpfr_mul(t1, a->im, b->re, ZRND);
    mpfr_mul(t2, a->re, b->im, ZRND);
    mpfr_sub(ni, t1, t2, ZRND);                  /* bc - ad */
    mpfr_div(out->re, nr, den, ZRND);
    mpfr_div(out->im, ni, den, ZRND);
    mpfr_clears(den, t1, t2, nr, ni, (mpfr_ptr)0);
}

static void zcx_exp(zcx* out, const zcx* z, mpfr_prec_t p) {
    mpfr_t ea, c, s;
    mpfr_inits2(p, ea, c, s, (mpfr_ptr)0);
    mpfr_exp(ea, z->re, ZRND);
    mpfr_sin_cos(s, c, z->im, ZRND);
    mpfr_mul(out->re, ea, c, ZRND);
    mpfr_mul(out->im, ea, s, ZRND);
    mpfr_clears(ea, c, s, (mpfr_ptr)0);
}

/* out = Log(z), principal branch. */
static void zcx_log(zcx* out, const zcx* z, mpfr_prec_t p) {
    mpfr_t mag, ang;
    mpfr_inits2(p, mag, ang, (mpfr_ptr)0);
    mpfr_hypot(mag, z->re, z->im, ZRND);
    mpfr_atan2(ang, z->im, z->re, ZRND);
    mpfr_log(out->re, mag, ZRND);
    mpfr_set(out->im, ang, ZRND);
    mpfr_clears(mag, ang, (mpfr_ptr)0);
}

/* out = base^expo (principal branch) = exp(expo * Log(base)). */
static void zcx_pow(zcx* out, const zcx* base, const zcx* expo, mpfr_prec_t p) {
    zcx lg, prod;
    zcx_init(&lg, p); zcx_init(&prod, p);
    zcx_log(&lg, base, p);
    zcx_mul(&prod, expo, &lg, p);
    zcx_exp(out, &prod, p);
    zcx_clear(&lg); zcx_clear(&prod);
}

/* ------------------------------------------------------------------ */
/* Euler-Maclaurin Hurwitz zeta kernel                                 */
/* ------------------------------------------------------------------ */

/* Number of head terms N for working precision wp and |s|. The correction
 * series error falls like (2 pi N)^-(2M); taking N ~ digits + |s| keeps M
 * modest and the head sum cheap. */
static unsigned long zeta_em_terms(mpfr_prec_t wp, const zcx* s) {
    double digits = (double)wp * 0.30103 + 5.0;
    double smag = hypot(mpfr_get_d(s->re, ZRND), mpfr_get_d(s->im, ZRND));
    double n = digits + smag + 10.0;
    if (n < 8.0) n = 8.0;
    if (n > 200000.0) n = 200000.0;
    return (unsigned long)n;
}

/* zeta(s, a) for complex s, a (Re(a) > 0 after the head shift) into `out`.
 * Returns false on a non-finite result (e.g. the s = 1 pole). */
static bool zeta_hurwitz_cx(zcx* out, const zcx* s, const zcx* a, mpfr_prec_t wp) {
    unsigned long N = zeta_em_terms(wp, s);

    zcx S, negs, base, term, w, P0, f, winv2, w2, sm1, r, tmp, tmp2;
    zcx_init(&S, wp); zcx_init(&negs, wp); zcx_init(&base, wp); zcx_init(&term, wp);
    zcx_init(&w, wp); zcx_init(&P0, wp); zcx_init(&f, wp); zcx_init(&winv2, wp);
    zcx_init(&w2, wp); zcx_init(&sm1, wp); zcx_init(&r, wp); zcx_init(&tmp, wp);
    zcx_init(&tmp2, wp);

    mpfr_t eps, mag, prevmag, sabs, cj, bq_d;
    mpfr_inits2(wp, eps, mag, prevmag, sabs, cj, bq_d, (mpfr_ptr)0);
    mpfr_set_ui(eps, 1, ZRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 12 ? wp - 8 : 1), ZRND);

    /* -s */
    mpfr_neg(negs.re, s->re, ZRND);
    mpfr_neg(negs.im, s->im, ZRND);

    /* Head sum S = Sum_{k=0}^{N-1} (a+k)^-s, skipping any k with a+k == 0. */
    mpfr_set_ui(S.re, 0, ZRND); mpfr_set_ui(S.im, 0, ZRND);
    for (unsigned long k = 0; k < N; k++) {
        mpfr_add_ui(base.re, a->re, k, ZRND);
        mpfr_set(base.im, a->im, ZRND);
        zcx_abs(mag, &base);
        if (mpfr_sgn(mag) == 0) continue;          /* k + a == 0 excluded */
        zcx_pow(&term, &base, &negs, wp);
        zcx_add(&S, &S, &term);
    }

    /* w = a + N. */
    mpfr_add_ui(w.re, a->re, N, ZRND);
    mpfr_set(w.im, a->im, ZRND);

    /* P0 = w^-s. */
    zcx_pow(&P0, &w, &negs, wp);

    /* + w^(1-s)/(s-1) = (P0 * w)/(s-1). */
    mpfr_sub_ui(sm1.re, s->re, 1, ZRND);
    mpfr_set(sm1.im, s->im, ZRND);
    zcx_mul(&tmp, &P0, &w, wp);
    zcx_div(&tmp2, &tmp, &sm1, wp);
    zcx_add(&S, &S, &tmp2);

    /* + 1/2 w^-s. */
    mpfr_set_d(cj, 0.5, ZRND);
    zcx_mul_r(&tmp, &P0, cj);
    zcx_add(&S, &S, &tmp);

    /* Correction series.  f_j = w^-(2j-1), starting f_1 = 1/w; advance by
     * winv2 = 1/w^2 each step.  r_j = (s)_{2j-1}, starting r_1 = s. */
    {
        zcx one; zcx_init(&one, wp);
        mpfr_set_ui(one.re, 1, ZRND); mpfr_set_ui(one.im, 0, ZRND);
        zcx_div(&f, &one, &w, wp);                 /* 1/w */
        zcx_mul(&w2, &w, &w, wp);
        zcx_div(&winv2, &one, &w2, wp);            /* 1/w^2 */
        zcx_clear(&one);
    }
    zcx_set(&r, s);                                /* (s)_1 = s */

    mpfr_set_inf(prevmag, 1);
    mpq_t bq;
    mpq_init(bq);
    mpz_t fac;
    mpz_init(fac);

    unsigned long jmax = (unsigned long)wp + 8;
    for (unsigned long j = 1; j <= jmax; j++) {
        /* coeff c_j = B_{2j}/(2j)!. */
        zbern_get_q(bq, (size_t)(2 * j));
        if (mpq_sgn(bq) == 0) {                    /* never (even index), but be safe */
            /* update r, f and continue */
        }
        zeta_factorial(fac, (long)(2 * j));
        mpfr_set_q(cj, bq, ZRND);
        mpfr_set_z(bq_d, fac, ZRND);
        mpfr_div(cj, cj, bq_d, ZRND);              /* B_{2j}/(2j)! */

        /* term = c_j * r * P0 * f. */
        zcx_mul(&tmp, &r, &P0, wp);
        zcx_mul(&tmp2, &tmp, &f, wp);
        zcx_mul_r(&term, &tmp2, cj);
        zcx_abs(mag, &term);

        /* Optimal truncation: stop once the terms start growing again. */
        if (j > 1 && mpfr_cmp(mag, prevmag) > 0) break;
        zcx_add(&S, &S, &term);
        mpfr_set(prevmag, mag, ZRND);

        /* Converged? |term| < eps * |S|. */
        zcx_abs(sabs, &S);
        mpfr_mul(sabs, sabs, eps, ZRND);
        if (mpfr_cmp(mag, sabs) < 0) break;

        /* r *= (s + 2j-1)(s + 2j). */
        mpfr_add_ui(tmp.re, s->re, 2 * j - 1, ZRND); mpfr_set(tmp.im, s->im, ZRND);
        zcx_mul(&tmp2, &r, &tmp, wp);
        mpfr_add_ui(tmp.re, s->re, 2 * j, ZRND);     mpfr_set(tmp.im, s->im, ZRND);
        zcx_mul(&r, &tmp2, &tmp, wp);
        /* f *= 1/w^2. */
        zcx_mul(&tmp, &f, &winv2, wp);
        zcx_set(&f, &tmp);
    }

    mpz_clear(fac);
    mpq_clear(bq);
    zcx_set(out, &S);
    bool ok = mpfr_number_p(out->re) && mpfr_number_p(out->im);

    mpfr_clears(eps, mag, prevmag, sabs, cj, bq_d, (mpfr_ptr)0);
    zcx_clear(&S); zcx_clear(&negs); zcx_clear(&base); zcx_clear(&term);
    zcx_clear(&w); zcx_clear(&P0); zcx_clear(&f); zcx_clear(&winv2);
    zcx_clear(&w2); zcx_clear(&sm1); zcx_clear(&r); zcx_clear(&tmp); zcx_clear(&tmp2);
    return ok;
}

/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool zeta_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, ZRND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          ZRND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        ZRND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          ZRND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, ZRND);
        mpfr_div_si(out, out, (long)d, ZRND);
        return true;
    }
    return false;
}

/* Result precision (bits) for Zeta[s, a] under Mathematica contagion: the
 * minimum precision among the inexact leaves of s and a, floored at machine
 * (53). A machine-precision argument forces a machine-precision result even if
 * the other argument is high-precision MPFR. See numeric_min_inexact_bits. */
static mpfr_prec_t zeta_out_prec(const Expr* s, const Expr* a) {
    long ps = numeric_min_inexact_bits(s);
    long pa = numeric_min_inexact_bits(a);
    long m  = (ps && pa) ? (ps < pa ? ps : pa) : (ps ? ps : pa);
    return m < 53 ? 53 : (mpfr_prec_t)m;
}

/* Build the numeric result: Real / Complex[Real,Real] for machine precision,
 * MPFR parts for arbitrary precision.  `real_only` forces a real result. */
static Expr* zeta_result(const mpfr_t re, const mpfr_t im,
                         mpfr_prec_t out_prec, bool real_only) {
    if (out_prec <= 53) {
        Expr* rr = expr_new_real(mpfr_get_d(re, ZRND));
        if (real_only) return rr;
        return make_complex(rr, expr_new_real(mpfr_get_d(im, ZRND)));
    }
    Expr* rr = expr_new_mpfr_bits(out_prec);
    mpfr_set(rr->data.mpfr, re, ZRND);
    if (real_only) return rr;
    Expr* ii = expr_new_mpfr_bits(out_prec);
    mpfr_set(ii->data.mpfr, im, ZRND);
    return make_complex(rr, ii);
}

/* Decompose a numeric value into borrowed (re, im). A plain real yields
 * (value, shared zero); a Complex[..] yields its two parts. */
static void zeta_decompose(Expr* v, Expr* zero, Expr** re, Expr** im) {
    Expr *r, *i;
    if (is_complex(v, &r, &i)) { *re = r; *im = i; }
    else { *re = v; *im = zero; }
}

/* Numeric Zeta[s, a] (a defaults to 1 for the Riemann case). At least one of
 * s, a must carry an inexact part. Returns Real/Complex/MPFR, ComplexInfinity
 * at the pole, or NULL if the parts are not numeric. */
static Expr* zeta_numeric(Expr* s, Expr* a) {
    Expr* zero = expr_new_integer(0);
    Expr *sre, *sim, *are, *aim;
    zeta_decompose(s, zero, &sre, &sim);
    zeta_decompose(a, zero, &are, &aim);

    mpfr_prec_t out_prec = zeta_out_prec(s, a);
    mpfr_prec_t wp = out_prec + 64;

    bool real_only = (sim == zero) && (aim == zero);

    zcx sc, ac, res;
    zcx_init(&sc, wp); zcx_init(&ac, wp); zcx_init(&res, wp);
    bool ok = zeta_set_mpfr(sc.re, sre) && zeta_set_mpfr(sc.im, sim) &&
              zeta_set_mpfr(ac.re, are) && zeta_set_mpfr(ac.im, aim);
    Expr* out = NULL;
    if (ok) {
        if (zeta_hurwitz_cx(&res, &sc, &ac, wp)) {
            out = zeta_result(res.re, res.im, out_prec, real_only);
        } else {
            out = expr_new_symbol("ComplexInfinity");
        }
    }
    zcx_clear(&sc); zcx_clear(&ac); zcx_clear(&res);
    expr_free(zero);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Zeta[s]  (Riemann zeta)                                             */
/* ------------------------------------------------------------------ */

static Expr* zeta_one_arg(Expr* s) {
    /* 1. Exact integer arguments -> closed forms. */
    long n;
    if (zeta_exact_int(s, &n)) {
        if (n == 1) return expr_new_symbol("ComplexInfinity");
        if (n == 0) return make_rational(-1, 2);
        if (n > 0) {
            if ((n & 1L) == 0) {                    /* even positive: rational*Pi^n */
                Expr* ev = zeta_exact_even(n);
                if (ev) return ev;
            }
            /* odd positive (>=3): no closed form -> symbolic (unless numeric) */
        } else {                                    /* n < 0 */
            Expr* nv = zeta_exact_negint(-n);
            if (nv) return nv;
        }
    }

    /* 2. Symbolic limiting value. */
    if (zeta_is_symbol(s, "Infinity")) return expr_new_integer(1);

#ifdef USE_MPFR
    /* 3/4. Numeric real Riemann zeta via MPFR (machine or arbitrary). */
    if (s->type == EXPR_REAL || s->type == EXPR_MPFR) {
        mpfr_prec_t prec = (s->type == EXPR_MPFR) ? mpfr_get_prec(s->data.mpfr) : 53;
        mpfr_t sv, rv;
        mpfr_init2(sv, prec + 16); mpfr_init2(rv, prec + 16);
        zeta_set_mpfr(sv, s);
        mpfr_zeta(rv, sv, ZRND);
        Expr* out;
        if (!mpfr_number_p(rv)) {
            out = expr_new_symbol("ComplexInfinity");   /* pole at s = 1 */
        } else if (s->type == EXPR_MPFR) {
            out = expr_new_mpfr_bits(prec);
            mpfr_set(out->data.mpfr, rv, ZRND);
        } else {
            out = expr_new_real(mpfr_get_d(rv, ZRND));
        }
        mpfr_clear(sv); mpfr_clear(rv);
        return out;
    }

    /* 5. Complex argument (Complex[..] with an inexact part) -> EM kernel
     *    with a = 1. */
    {
        Expr *re, *im;
        if (is_complex(s, &re, &im) &&
            (zeta_is_inexact(re) || zeta_is_inexact(im))) {
            Expr* one = expr_new_integer(1);
            Expr* out = zeta_numeric(s, one);
            expr_free(one);
            if (out) return out;
        }
    }
#endif

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Zeta[s, a]  (Hurwitz zeta)                                          */
/* ------------------------------------------------------------------ */

/* Exact Hurwitz at a positive integer a = m:
 *   zeta(s, m) = Zeta[s] - Sum_{k=1}^{m-1} k^-s.
 * Valid for symbolic / exact / integer s (gives Zeta[3,2] -> -1+Zeta[3],
 * Zeta[4,5] -> Pi^4/90 - 22369/20736). Returns NULL past the cap. */
static Expr* zeta_hurwitz_int_a(Expr* s, long m) {
    if (m < 1 || m > ZETA_HURWITZ_A_CAP) return NULL;

    Expr* zs = expr_new_function(expr_new_symbol("Zeta"), (Expr*[]){ expr_copy(s) }, 1);
    if (m == 1) return eval_and_free(zs);          /* zeta(s, 1) = zeta(s) */

    /* Build Zeta[s] - Sum_{k=1}^{m-1} k^-s. */
    size_t terms = (size_t)m;                       /* zs + (m-1) subtracted terms */
    Expr** parts = (Expr**)malloc(terms * sizeof(Expr*));
    if (!parts) { expr_free(zs); return NULL; }
    parts[0] = zs;
    for (long k = 1; k < m; k++) {
        /* -(k^-s) = Times[-1, Power[k, -s]]. */
        Expr* negs = expr_new_function(expr_new_symbol("Times"),
                         (Expr*[]){ expr_new_integer(-1), expr_copy(s) }, 2);
        Expr* kpow = expr_new_function(expr_new_symbol("Power"),
                         (Expr*[]){ expr_new_integer(k), negs }, 2);
        parts[(size_t)k] = expr_new_function(expr_new_symbol("Times"),
                         (Expr*[]){ expr_new_integer(-1), kpow }, 2);
    }
    Expr* sum = expr_new_function(expr_new_symbol("Plus"), parts, terms);
    free(parts);
    return eval_and_free(sum);
}

static Expr* zeta_two_arg(Expr* s, Expr* a) {
    /* 1. zeta(s, 1) = zeta(s). */
    if (a->type == EXPR_INTEGER && a->data.integer == 1)
        return eval_and_free(expr_new_function(expr_new_symbol("Zeta"),
                             (Expr*[]){ expr_copy(s) }, 1));

    /* 2. Exact Hurwitz at a positive integer a (only when s is not inexact;
     *    inexact s falls through to the numeric kernel). */
    if (!zeta_is_inexact(s)) {
        long m;
        if (zeta_exact_int(a, &m) && m >= 1) {
            Expr* closed = zeta_hurwitz_int_a(s, m);
            if (closed) return closed;
        }
    }

#ifdef USE_MPFR
    /* 3. Numeric (real or complex) Hurwitz: at least one inexact operand. */
    {
        Expr* zero = expr_new_integer(0);
        Expr *sre, *sim, *are, *aim;
        zeta_decompose(s, zero, &sre, &sim);
        zeta_decompose(a, zero, &are, &aim);
        bool any_inexact = zeta_is_inexact(sre) || zeta_is_inexact(sim) ||
                           zeta_is_inexact(are) || zeta_is_inexact(aim);
        expr_free(zero);
        if (any_inexact) {
            Expr* out = zeta_numeric(s, a);
            if (out) return out;
        }
    }
#endif

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

static Expr* zeta_emit_argt(size_t argc) {
    fprintf(stderr,
            "Zeta::argt: Zeta called with %zu argument%s; "
            "1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_zeta(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return zeta_one_arg(args[0]);
    if (argc == 2) return zeta_two_arg(args[0], args[1]);
    return zeta_emit_argt(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void zeta_init(void) {
    symtab_add_builtin("Zeta", builtin_zeta);
    symtab_get_def("Zeta")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

/* Mathilda -- the PolyGamma function family.
 *
 *   PolyGamma[z]      digamma psi(z) = Gamma'(z)/Gamma(z)
 *                       (always rewrites to the two-argument form PolyGamma[0, z])
 *   PolyGamma[n, z]   n-th polygamma psi^(n)(z) = d^n/dz^n psi(z)
 *
 * Evaluation is layered so each kind of argument takes the cheapest exact
 * or fastest numeric route available:
 *
 *   order n = -1                -> LogGamma[z]   (inert; psi^(-1) = log-gamma)
 *   z a non-positive integer    -> ComplexInfinity (pole of every psi^(n))
 *   z a positive integer        -> exact closed form:
 *                                    n = 0          : H_{z-1} - EulerGamma
 *                                    n >= 1 odd     : rational + rational*Pi^(n+1)
 *                                                     (via zeta(n+1), even -> Bernoulli)
 *                                    n >= 1 even    : stays symbolic (zeta(odd), no
 *                                                     closed form)
 *   z inexact real (Real/MPFR)  -> numeric: mpfr_digamma for n = 0, otherwise the
 *                                    recurrence-shift + Bernoulli asymptotic series
 *   z inexact complex           -> the same asymptotic, in complex arithmetic
 *   everything else             -> stays symbolic (return NULL)
 *
 * MPFR provides a digamma (n = 0) but no higher polygamma and no Hurwitz zeta,
 * so the n >= 1 numeric kernel and the complex digamma are implemented here from
 * the classical asymptotic expansion (Abramowitz & Stegun 6.4.11):
 *
 *   psi^(n)(w) ~ (-1)^(n-1) [ (n-1)!/w^n + n!/(2 w^(n+1))
 *                             + Sum_{k>=1} B_{2k} (2k+n-1)!/(2k)! w^(-(2k+n)) ]   (n>=1)
 *   psi(w)     ~ ln w - 1/(2w) - Sum_{k>=1} B_{2k}/(2k) w^(-2k)                   (n=0)
 *
 * valid for large |w|; the recurrence psi^(n)(z) = psi^(n)(z+1) - (-1)^n n! z^(-(n+1))
 * shifts a small or negative argument up into that regime.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "polygamma.h"
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
#define GRND MPFR_RNDN
#endif

/* Largest positive integer argument for which we build the exact closed form;
 * beyond this the harmonic / power sums get unwieldy, so stay symbolic. */
#define POLYGAMMA_EXACT_ARG_CAP   100000L
/* Largest order for which we attempt the exact (Bernoulli) closed form. */
#define POLYGAMMA_EXACT_ORDER_CAP 256L
/* Largest order for which we attempt numeric evaluation at all. */
#define POLYGAMMA_NUMERIC_ORDER_CAP 1024L

/* ------------------------------------------------------------------ */
/* Small predicates / coercions                                        */
/* ------------------------------------------------------------------ */

/* True if `e` is exactly the symbol `name`. */
static bool pg_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol, name) == 0;
}

#ifdef USE_MPFR
/* True if `e` is an inexact numeric leaf (Real or MPFR). Only the numeric
 * (MPFR-gated) dispatch path consults this. */
static bool pg_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
    if (e->type == EXPR_MPFR) return true;
    return false;
}
#endif

/* Extract a machine integer order. Accepts exact Integer / BigInt, and also
 * integer-valued inexact reals (Real / MPFR) -- N[] numericalises the order
 * argument (0 -> 0.0), and PolyGamma[0., z] must still mean the digamma. */
static bool pg_small_int(const Expr* e, long* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (long)e->data.integer; return true; }
    if (e->type == EXPR_BIGINT && mpz_fits_slong_p(e->data.bigint)) {
        *out = mpz_get_si(e->data.bigint);
        return true;
    }
    if (e->type == EXPR_REAL) {
        double v = e->data.real;
        if (v == floor(v) && fabs(v) < 1e15) { *out = (long)v; return true; }
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR && mpfr_integer_p(e->data.mpfr) &&
        mpfr_fits_slong_p(e->data.mpfr, GRND)) {
        *out = mpfr_get_si(e->data.mpfr, GRND);
        return true;
    }
#endif
    return false;
}

/* True if `e` is a non-positive exact integer (a pole of every psi^(n)). */
static bool pg_nonpos_int(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer <= 0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) <= 0;
    return false;
}

/* True (and sets *m) if `e` is a positive exact integer fitting a long. */
static bool pg_pos_int(const Expr* e, long* m) {
    if (e->type == EXPR_INTEGER && e->data.integer > 0) {
        *m = (long)e->data.integer;
        return true;
    }
    if (e->type == EXPR_BIGINT && mpz_sgn(e->data.bigint) > 0 &&
        mpz_fits_slong_p(e->data.bigint)) {
        *m = mpz_get_si(e->data.bigint);
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ */
/* Exact rational -> Expr helpers                                      */
/* ------------------------------------------------------------------ */

/* mpz -> Integer/BigInt Expr (demoted to int64 when it fits). */
static Expr* expr_from_mpz(const mpz_t z) {
    return expr_bigint_normalize(expr_new_bigint_from_mpz(z));
}

/* mpq -> exact Expr: an Integer/BigInt for whole values, else the canonical
 * Rational[num, den] (Mathilda's internal representation of a fraction). */
static Expr* expr_from_mpq(const mpq_t q) {
    mpq_t t;
    mpq_init(t);
    mpq_set(t, q);
    mpq_canonicalize(t);
    Expr* out;
    if (mpz_cmp_ui(mpq_denref(t), 1) == 0) {
        out = expr_from_mpz(mpq_numref(t));
    } else {
        Expr* num = expr_from_mpz(mpq_numref(t));
        Expr* den = expr_from_mpz(mpq_denref(t));
        Expr* args[2] = { num, den };
        out = expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
    }
    mpq_clear(t);
    return out;
}

/* ------------------------------------------------------------------ */
/* Bernoulli numbers B_k (exact, lazily cached)                        */
/*                                                                     */
/* Recurrence: B_0 = 1, and for m >= 1                                 */
/*   B_m = -1/(m+1) Sum_{k=0}^{m-1} C(m+1, k) B_k.                      */
/* Only even-index Bernoulli numbers (and B_0, B_1) are ever consumed. */
/* The cache is a process-lifetime singleton (reachable, not leaked).  */
/* ------------------------------------------------------------------ */

static mpq_t* g_bern = NULL;
static size_t g_bern_len = 0;

static void bern_ensure(size_t upto) {
    if (upto + 1 <= g_bern_len) return;
    size_t newlen = upto + 1;
    mpq_t* grown = (mpq_t*)realloc(g_bern, newlen * sizeof(mpq_t));
    if (!grown) return; /* out of memory: leave cache as-is, callers degrade */
    g_bern = grown;

    mpq_t sum, term, factor;
    mpq_inits(sum, term, factor, (mpq_ptr)0);
    mpz_t binz;
    mpz_init(binz);

    for (size_t m = g_bern_len; m < newlen; m++) {
        mpq_init(g_bern[m]);
        if (m == 0) {
            mpq_set_ui(g_bern[m], 1, 1);
            continue;
        }
        mpq_set_ui(sum, 0, 1);
        for (size_t k = 0; k < m; k++) {
            mpz_bin_uiui(binz, (unsigned long)(m + 1), (unsigned long)k);
            mpq_set_z(factor, binz);          /* C(m+1, k) */
            mpq_mul(term, factor, g_bern[k]);
            mpq_add(sum, sum, term);
        }
        mpq_set_si(term, -1, (unsigned long)(m + 1)); /* -1/(m+1) */
        mpq_canonicalize(term);
        mpq_mul(g_bern[m], sum, term);
        mpq_canonicalize(g_bern[m]);
    }
    g_bern_len = newlen;

    mpz_clear(binz);
    mpq_clears(sum, term, factor, (mpq_ptr)0);
}

/* B_idx into `out`. Odd indices above 1 are exactly zero. */
static void bern_get_q(mpq_t out, size_t idx) {
    if (idx > 1 && (idx & 1u)) { mpq_set_ui(out, 0, 1); return; }
    bern_ensure(idx);
    if (idx < g_bern_len) mpq_set(out, g_bern[idx]);
    else mpq_set_ui(out, 0, 1); /* allocation failure fallback */
}

/* ------------------------------------------------------------------ */
/* Exact closed forms at positive integer arguments                    */
/* ------------------------------------------------------------------ */

/* H_{m-1} = Sum_{k=1}^{m-1} 1/k. */
static void harmonic_sum(mpq_t out, long m) {
    mpq_set_ui(out, 0, 1);
    mpq_t t;
    mpq_init(t);
    for (long k = 1; k < m; k++) {
        mpq_set_ui(t, 1, (unsigned long)k);
        mpq_add(out, out, t);
    }
    mpq_clear(t);
}

/* Sum_{k=1}^{m-1} 1/k^p  (p >= 1). */
static void power_sum(mpq_t out, long m, long p) {
    mpq_set_ui(out, 0, 1);
    mpq_t t;
    mpz_t kp;
    mpq_init(t);
    mpz_init(kp);
    for (long k = 1; k < m; k++) {
        mpz_ui_pow_ui(kp, (unsigned long)k, (unsigned long)p);
        mpq_set_z(t, kp);          /* k^p */
        mpq_inv(t, t);             /* 1/k^p */
        mpq_add(out, out, t);
    }
    mpz_clear(kp);
    mpq_clear(t);
}

static void fact_mpz(mpz_t out, long n) {
    mpz_set_ui(out, 1);
    for (long i = 2; i <= n; i++) mpz_mul_ui(out, out, (unsigned long)i);
}

/* Exact PolyGamma[n, m] for integer order n >= 0 and positive integer m.
 * Returns an owned Expr, or NULL when there is no closed form we emit
 * (even order n >= 2, or arguments past the caps). */
static Expr* pg_exact_at_positive_int(long n, long m) {
    if (m > POLYGAMMA_EXACT_ARG_CAP) return NULL;

    if (n == 0) {
        /* psi(m) = H_{m-1} - EulerGamma. */
        mpq_t H;
        mpq_init(H);
        harmonic_sum(H, m);
        Expr* hrat = expr_from_mpq(H);
        mpq_clear(H);
        Expr* eg_args[2] = { expr_new_integer(-1), expr_new_symbol(SYM_EulerGamma) };
        Expr* eg = expr_new_function(expr_new_symbol(SYM_Times), eg_args, 2);
        Expr* sum_args[2] = { hrat, eg };
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), sum_args, 2);
        return eval_and_free(sum);
    }

    if ((n % 2) == 0) return NULL;           /* even order: zeta(odd), no closed form */
    if (n > POLYGAMMA_EXACT_ORDER_CAP) return NULL;

    /* psi^(n)(m) = (-1)^(n+1) n! [ zeta(n+1) - Sum_{k=1}^{m-1} 1/k^(n+1) ].
     * For odd n, n+1 is even so (-1)^(n+1) = +1 and
     *   zeta(n+1) = c * Pi^(n+1),  c = (-1)^(j+1) B_{2j} 2^(2j-1) / (2j)!,  2j = n+1.
     * Hence the result is  A*Pi^(n+1) + B  with
     *   A = n! c,   B = -n! Sum 1/k^(n+1). */
    long j = (n + 1) / 2;

    mpq_t c, tmpq;
    mpz_t pw, f2j, nf;
    mpq_inits(c, tmpq, (mpq_ptr)0);
    mpz_inits(pw, f2j, nf, (mpz_ptr)0);

    bern_get_q(c, (size_t)(2 * j));                 /* B_{2j} */
    mpz_ui_pow_ui(pw, 2, (unsigned long)(2 * j - 1)); /* 2^(2j-1) */
    mpq_set_z(tmpq, pw);
    mpq_mul(c, c, tmpq);
    fact_mpz(f2j, 2 * j);                            /* (2j)! */
    mpq_set_z(tmpq, f2j);
    mpq_div(c, c, tmpq);
    if ((j % 2) == 0) mpq_neg(c, c);                 /* (-1)^(j+1) */

    fact_mpz(nf, n);                                 /* n! */
    mpq_set_z(tmpq, nf);

    mpq_t A, B, S;
    mpq_inits(A, B, S, (mpq_ptr)0);
    mpq_mul(A, tmpq, c);                             /* A = n! c */
    mpq_canonicalize(A);
    power_sum(S, m, n + 1);
    mpq_mul(B, tmpq, S);                             /* n! S */
    mpq_neg(B, B);                                   /* B = -n! S */
    mpq_canonicalize(B);

    Expr* Aex = expr_from_mpq(A);
    Expr* pipow_args[2] = { expr_new_symbol(SYM_Pi), expr_new_integer(n + 1) };
    Expr* pipow = expr_new_function(expr_new_symbol(SYM_Power), pipow_args, 2);
    Expr* term_args[2] = { Aex, pipow };
    Expr* term = expr_new_function(expr_new_symbol(SYM_Times), term_args, 2);
    Expr* Bex = expr_from_mpq(B);
    Expr* sum_args[2] = { term, Bex };
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), sum_args, 2);

    mpq_clears(c, tmpq, A, B, S, (mpq_ptr)0);
    mpz_clears(pw, f2j, nf, (mpz_ptr)0);

    return eval_and_free(sum);
}

/* ================================================================== */
/* Numeric evaluation (MPFR)                                           */
/* ================================================================== */
#ifdef USE_MPFR

/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool pg_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, GRND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          GRND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        GRND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          GRND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, GRND);
        mpfr_div_si(out, out, (long)d, GRND);
        return true;
    }
    return false;
}

/* ---- minimal complex-MPFR toolkit (re/im pair); ops are alias-safe ---- */
typedef struct { mpfr_t re, im; } pcx;

static void pcx_init(pcx* z, mpfr_prec_t p) { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void pcx_clear(pcx* z)               { mpfr_clear(z->re);    mpfr_clear(z->im);    }
static void pcx_set(pcx* d, const pcx* s)   { mpfr_set(d->re, s->re, GRND); mpfr_set(d->im, s->im, GRND); }
static void pcx_add(pcx* o, const pcx* a, const pcx* b) {
    mpfr_add(o->re, a->re, b->re, GRND);
    mpfr_add(o->im, a->im, b->im, GRND);
}
static void pcx_sub(pcx* o, const pcx* a, const pcx* b) {
    mpfr_sub(o->re, a->re, b->re, GRND);
    mpfr_sub(o->im, a->im, b->im, GRND);
}
static void pcx_abs(mpfr_t mag, const pcx* z) { mpfr_hypot(mag, z->re, z->im, GRND); }

static void pcx_mul(pcx* o, const pcx* a, const pcx* b, mpfr_prec_t p) {
    mpfr_t ac, bd, ad, bc;
    mpfr_inits2(p, ac, bd, ad, bc, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, GRND);
    mpfr_mul(bd, a->im, b->im, GRND);
    mpfr_mul(ad, a->re, b->im, GRND);
    mpfr_mul(bc, a->im, b->re, GRND);
    mpfr_sub(o->re, ac, bd, GRND);
    mpfr_add(o->im, ad, bc, GRND);
    mpfr_clears(ac, bd, ad, bc, (mpfr_ptr)0);
}

static void pcx_div(pcx* o, const pcx* a, const pcx* b, mpfr_prec_t p) {
    mpfr_t den, t1, t2, nr, ni;
    mpfr_inits2(p, den, t1, t2, nr, ni, (mpfr_ptr)0);
    mpfr_mul(t1, b->re, b->re, GRND);
    mpfr_mul(t2, b->im, b->im, GRND);
    mpfr_add(den, t1, t2, GRND);
    mpfr_mul(t1, a->re, b->re, GRND);
    mpfr_mul(t2, a->im, b->im, GRND);
    mpfr_add(nr, t1, t2, GRND);          /* ac + bd */
    mpfr_mul(t1, a->im, b->re, GRND);
    mpfr_mul(t2, a->re, b->im, GRND);
    mpfr_sub(ni, t1, t2, GRND);          /* bc - ad */
    mpfr_div(o->re, nr, den, GRND);
    mpfr_div(o->im, ni, den, GRND);
    mpfr_clears(den, t1, t2, nr, ni, (mpfr_ptr)0);
}

/* o = a * s  (real scalar). */
static void pcx_mul_fr(pcx* o, const pcx* a, const mpfr_t s) {
    mpfr_mul(o->re, a->re, s, GRND);
    mpfr_mul(o->im, a->im, s, GRND);
}

/* o = Log(a), principal branch. */
static void pcx_log(pcx* o, const pcx* a, mpfr_prec_t p) {
    mpfr_t mag, ang;
    mpfr_inits2(p, mag, ang, (mpfr_ptr)0);
    mpfr_hypot(mag, a->re, a->im, GRND);
    mpfr_atan2(ang, a->im, a->re, GRND);
    mpfr_log(o->re, mag, GRND);
    mpfr_set(o->im, ang, GRND);
    mpfr_clears(mag, ang, (mpfr_ptr)0);
}

/* Core: psi^(n)(z) for integer n >= 0 and complex z, to `prec` bits.
 * Writes the result into (outre, outim); returns false if it could not
 * compute a finite value (NaN/Inf, or a pathological shift). */
static bool pcx_polygamma(long n, const mpfr_t zre, const mpfr_t zim,
                          mpfr_t outre, mpfr_t outim, mpfr_prec_t prec) {
    if (!mpfr_number_p(zre) || !mpfr_number_p(zim)) return false;

    mpfr_prec_t wp = prec + 64 + (mpfr_prec_t)(2 * (n + 2));

    /* Shift count M so the asymptotic series converges: push Re(z) up past
     * ~0.13*wp + n.  The smallest-term test then truncates the series. */
    double zr = mpfr_get_d(zre, GRND);
    long target = (long)(0.13 * (double)wp) + n + 8;
    long M = target - (long)floor(zr);
    if (M < 0) M = 0;
    if (M > 2000000L) return false;

    pcx z, w, r, r2, one, sum, pk, term;
    pcx_init(&z, wp);  pcx_init(&w, wp);  pcx_init(&r, wp);   pcx_init(&r2, wp);
    pcx_init(&one, wp);pcx_init(&sum, wp);pcx_init(&pk, wp);  pcx_init(&term, wp);
    mpfr_t eps, mag, prevmag, coeff;
    mpfr_inits2(wp, eps, mag, prevmag, coeff, (mpfr_ptr)0);

    mpfr_set(z.re, zre, GRND); mpfr_set(z.im, zim, GRND);
    pcx_set(&w, &z); mpfr_add_si(w.re, w.re, M, GRND);   /* w = z + M */
    mpfr_set_si(one.re, 1, GRND); mpfr_set_zero(one.im, 1);
    pcx_div(&r, &one, &w, wp);                            /* r  = 1/w     */
    pcx_mul(&r2, &r, &r, wp);                             /* r2 = 1/w^2   */

    mpfr_set_ui(eps, 1, GRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(prec + 8), GRND);

    mpq_t bq;
    mpq_init(bq);

    if (n == 0) {
        /* psi(w) = log(w) - 1/(2w) - Sum_{k>=1} B_{2k}/(2k) w^(-2k). */
        pcx_log(&sum, &w, wp);
        mpfr_set_d(coeff, 0.5, GRND);
        pcx_mul_fr(&term, &r, coeff);
        pcx_sub(&sum, &sum, &term);
        pcx_set(&pk, &r2);                               /* w^(-2k), k=1 */
        mpfr_set_inf(prevmag, 1);
        for (long k = 1; k <= (long)wp; k++) {
            bern_get_q(bq, (size_t)(2 * k));
            mpfr_set_q(coeff, bq, GRND);
            mpfr_div_ui(coeff, coeff, (unsigned long)(2 * k), GRND);
            pcx_mul_fr(&term, &pk, coeff);
            pcx_abs(mag, &term);
            if (mpfr_cmp(mag, prevmag) > 0) break;        /* past optimal truncation */
            pcx_sub(&sum, &sum, &term);
            mpfr_set(prevmag, mag, GRND);
            if (mpfr_cmp(mag, eps) < 0) break;
            pcx_mul(&pk, &pk, &r2, wp);
        }
    } else {
        /* psi^(n)(w) = (-1)^(n-1) [ (n-1)!/w^n + n!/(2 w^(n+1))
         *               + Sum_{k>=1} B_{2k} (2k+n-1)!/(2k)! w^(-(2k+n)) ]. */
        int neg = (n % 2 == 0);                           /* (-1)^(n-1) = -1 iff n even */
        pcx rn, rn1, tmp;
        pcx_init(&rn, wp); pcx_init(&rn1, wp); pcx_init(&tmp, wp);
        pcx_set(&rn, &r);
        for (long i = 1; i < n; i++) pcx_mul(&rn, &rn, &r, wp); /* rn = w^(-n) */
        pcx_mul(&rn1, &rn, &r, wp);                            /* rn1 = w^(-(n+1)) */

        mpz_t f;
        mpz_init(f);
        fact_mpz(f, n - 1);                              /* (n-1)! */
        mpfr_set_z(coeff, f, GRND);
        pcx_mul_fr(&sum, &rn, coeff);                    /* (n-1)! w^(-n) */
        mpz_mul_ui(f, f, (unsigned long)n);              /* n! */
        mpfr_set_z(coeff, f, GRND);
        mpfr_div_ui(coeff, coeff, 2, GRND);
        pcx_mul_fr(&tmp, &rn1, coeff);                   /* n!/2 w^(-(n+1)) */
        pcx_add(&sum, &sum, &tmp);

        pcx_set(&pk, &r2);                               /* (1/w^2)^k, k=1 */
        mpfr_set_inf(prevmag, 1);
        for (long k = 1; k <= (long)wp; k++) {
            bern_get_q(bq, (size_t)(2 * k));
            mpfr_set_q(coeff, bq, GRND);                 /* B_{2k} */
            for (long i = 2 * k + 1; i <= 2 * k + n - 1; i++)
                mpfr_mul_ui(coeff, coeff, (unsigned long)i, GRND); /* *(2k+n-1)!/(2k)! */
            pcx_mul(&term, &rn, &pk, wp);                /* w^(-(2k+n)) */
            pcx_mul_fr(&term, &term, coeff);
            pcx_abs(mag, &term);
            if (mpfr_cmp(mag, prevmag) > 0) break;
            pcx_add(&sum, &sum, &term);
            mpfr_set(prevmag, mag, GRND);
            if (mpfr_cmp(mag, eps) < 0) break;
            pcx_mul(&pk, &pk, &r2, wp);
        }
        if (neg) { mpfr_neg(sum.re, sum.re, GRND); mpfr_neg(sum.im, sum.im, GRND); }

        mpz_clear(f);
        pcx_clear(&rn); pcx_clear(&rn1); pcx_clear(&tmp);
    }

    /* Undo the shift:
     *   psi^(n)(z) = psi^(n)(w) - (-1)^n n! Sum_{j=0}^{M-1} (z+j)^(-(n+1)). */
    if (M > 0) {
        pcx corr, u, up, base;
        pcx_init(&corr, wp); pcx_init(&u, wp); pcx_init(&up, wp); pcx_init(&base, wp);
        mpfr_set_zero(corr.re, 1); mpfr_set_zero(corr.im, 1);
        for (long jj = 0; jj < M; jj++) {
            pcx_set(&u, &z); mpfr_add_si(u.re, u.re, jj, GRND);
            pcx_div(&up, &one, &u, wp);                  /* 1/(z+j) */
            pcx_set(&base, &up);
            for (long i = 1; i < n + 1; i++) pcx_mul(&up, &up, &base, wp); /* (z+j)^(-(n+1)) */
            pcx_add(&corr, &corr, &up);
        }
        mpz_t nf;
        mpz_init(nf);
        fact_mpz(nf, n);                                 /* n! */
        mpfr_set_z(coeff, nf, GRND);
        if (n % 2 == 1) mpfr_neg(coeff, coeff, GRND);    /* (-1)^n */
        mpz_clear(nf);
        pcx_mul_fr(&corr, &corr, coeff);
        pcx_sub(&sum, &sum, &corr);
        pcx_clear(&corr); pcx_clear(&u); pcx_clear(&up); pcx_clear(&base);
    }

    bool ok = mpfr_number_p(sum.re) && mpfr_number_p(sum.im);
    if (ok) { mpfr_set(outre, sum.re, GRND); mpfr_set(outim, sum.im, GRND); }

    mpq_clear(bq);
    pcx_clear(&z); pcx_clear(&w); pcx_clear(&r); pcx_clear(&r2);
    pcx_clear(&one); pcx_clear(&sum); pcx_clear(&pk); pcx_clear(&term);
    mpfr_clears(eps, mag, prevmag, coeff, (mpfr_ptr)0);
    return ok;
}

/* PolyGamma[n, z] for an inexact real z (Real or MPFR). */
static Expr* pg_numeric_real(long n, const Expr* z) {
    mpfr_prec_t prec = (z->type == EXPR_MPFR) ? mpfr_get_prec(z->data.mpfr) : 53;
    mpfr_t x;
    mpfr_init2(x, prec);
    if (!pg_set_mpfr(x, z)) { mpfr_clear(x); return NULL; }

    Expr* out = NULL;
    if (n == 0) {
        mpfr_t rv;
        mpfr_init2(rv, prec);
        mpfr_digamma(rv, x, GRND);
        if (mpfr_nan_p(rv))      out = NULL;
        else if (mpfr_inf_p(rv)) out = expr_new_symbol(SYM_ComplexInfinity);
        else out = (z->type == EXPR_MPFR) ? expr_new_mpfr_copy(rv)
                                          : expr_new_real(mpfr_get_d(rv, GRND));
        mpfr_clear(rv);
    } else {
        mpfr_t zi, rre, rim;
        mpfr_inits2(prec, zi, rre, rim, (mpfr_ptr)0);
        mpfr_set_zero(zi, 1);
        if (pcx_polygamma(n, x, zi, rre, rim, prec))
            out = (z->type == EXPR_MPFR) ? expr_new_mpfr_copy(rre)
                                         : expr_new_real(mpfr_get_d(rre, GRND));
        mpfr_clears(zi, rre, rim, (mpfr_ptr)0);
    }
    mpfr_clear(x);
    return out;
}

/* PolyGamma[n, Complex[re, im]] with at least one inexact part. */
static Expr* pg_numeric_complex(long n, Expr* re, Expr* im) {
    /* Min-contagion across the inexact parts, floored at machine (53). A
     * machine-precision part forces a machine-precision result. */
    long pr = numeric_min_inexact_bits(re);
    long pi = numeric_min_inexact_bits(im);
    long mbits = (pr && pi) ? (pr < pi ? pr : pi) : (pr ? pr : pi);
    mpfr_prec_t prec = mbits < 53 ? 53 : (mpfr_prec_t)mbits;
    bool any_mpfr = (prec > 53);

    mpfr_t xr, xi, rr, ri;
    mpfr_inits2(prec, xr, xi, rr, ri, (mpfr_ptr)0);
    Expr* out = NULL;
    if (pg_set_mpfr(xr, re) && pg_set_mpfr(xi, im) &&
        pcx_polygamma(n, xr, xi, rr, ri, prec)) {
        Expr *er, *ei;
        if (any_mpfr) { er = expr_new_mpfr_copy(rr); ei = expr_new_mpfr_copy(ri); }
        else          { er = expr_new_real(mpfr_get_d(rr, GRND));
                        ei = expr_new_real(mpfr_get_d(ri, GRND)); }
        out = make_complex(er, ei);
    }
    mpfr_clears(xr, xi, rr, ri, (mpfr_ptr)0);
    return out;
}
#endif /* USE_MPFR */

/* ================================================================== */
/* Dispatch                                                            */
/* ================================================================== */

static Expr* polygamma_two_arg(Expr* order, Expr* z) {
    long n;
    if (!pg_small_int(order, &n)) return NULL; /* symbolic / fractional / complex order */

    if (n == -1) {
        /* psi^(-1) = log-gamma; emitted as the inert LogGamma. */
        Expr* zc = expr_copy(z);
        return expr_new_function(expr_new_symbol(SYM_LogGamma), &zc, 1);
    }
    if (n < 0) return NULL;                     /* orders <= -2 stay symbolic */

    /* Special values. */
    if (pg_is_symbol(z, "ComplexInfinity")) return expr_new_symbol(SYM_Indeterminate);
    if (pg_is_symbol(z, "Indeterminate"))   return expr_new_symbol(SYM_Indeterminate);
    if (pg_is_symbol(z, "Infinity"))
        return (n == 0) ? expr_new_symbol(SYM_Infinity) : expr_new_integer(0);
    if (pg_nonpos_int(z)) return expr_new_symbol(SYM_ComplexInfinity); /* pole */

    /* Exact closed form at positive integer arguments. */
    long m;
    if (pg_pos_int(z, &m)) {
        Expr* exact = pg_exact_at_positive_int(n, m);
        return exact; /* NULL -> stays symbolic (even order, or past caps) */
    }

    if (n > POLYGAMMA_NUMERIC_ORDER_CAP) return NULL;

#ifdef USE_MPFR
    /* Numeric: inexact real, then inexact complex. */
    if (pg_is_inexact(z)) return pg_numeric_real(n, z);
    Expr *re, *im;
    if (is_complex(z, &re, &im) && (pg_is_inexact(re) || pg_is_inexact(im)))
        return pg_numeric_complex(n, re, im);
#endif

    return NULL; /* leave symbolic */
}

static Expr* pg_emit_argt(size_t argc) {
    fprintf(stderr,
            "PolyGamma::argt: PolyGamma called with %zu argument%s; "
            "1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_polygamma(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) {
        /* PolyGamma[z] -> PolyGamma[0, z] (re-evaluated to a fixed point). */
        Expr* pg_args[2] = { expr_new_integer(0), expr_copy(args[0]) };
        return expr_new_function(expr_new_symbol(SYM_PolyGamma), pg_args, 2);
    }
    if (argc == 2) return polygamma_two_arg(args[0], args[1]);
    return pg_emit_argt(argc);
}

void polygamma_init(void) {
    symtab_add_builtin("PolyGamma", builtin_polygamma);
    symtab_get_def("PolyGamma")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);

    /* PolyGamma[-1, z] emits LogGamma[z]; the LogGamma builtin itself is
     * registered separately by loggamma_init (see loggamma.c). */
}

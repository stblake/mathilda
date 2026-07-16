/* Mathilda -- LogGamma[z], the log-gamma function log(Gamma(z)).
 *
 * LogGamma is analytic throughout the complex plane except for a branch cut on
 * the negative real axis. It is the analytic continuation of log(Gamma(z)) and
 * is *not* Log[Gamma[z]] (which has a more involved branch structure). Because
 * it never overflows where Gamma does, it is the natural primitive for
 * factorial ratios and asymptotics.
 *
 * Evaluation is layered, mirroring Gamma (gamma.c), so each kind of argument
 * takes the cheapest exact or fastest numeric route:
 *
 *   exact integer n >= 1          ->  Log[(n-1)!]   (exact, via Gamma's machinery)
 *   exact half-integer            ->  Log of the exact Sqrt[Pi] form, with the
 *                                     branch term -Ceiling[-z] Pi I for z < 0
 *   non-positive integer          ->  Infinity      (pole)
 *   symbolic infinities           ->  Infinity / Indeterminate / ComplexInfinity
 *   machine real z > 0            ->  libm   lgamma
 *   machine real z < 0 (non-int)  ->  lgamma + branch term (complex result)
 *   arbitrary real                ->  MPFR   mpfr_lgamma  (+ branch term)
 *   machine complex               ->  Lanczos log-gamma (double complex)
 *   arbitrary complex             ->  Spouge log-gamma   (MPFR, runtime coeffs)
 *   everything else               ->  stays symbolic (return NULL)
 *
 * The imaginary branch term on the negative real axis is Im = -Pi Ceiling[-z],
 * taken from above (Mathematica's convention); this is exact for the symbolic
 * half-integer path and is reused by the numeric real paths.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "loggamma.h"
#include "sym_names.h"

#include <complex.h>
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
#endif

/* M_PI is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Small numeric / structural helpers                                 */
/* ------------------------------------------------------------------ */

/* Coerce an exact-or-real leaf to a double. Succeeds for Integer, Real,
 * BigInt and Rational; fails for symbols, complex, MPFR and anything else. */
static bool lg_to_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { return false; } /* handled by the MPFR path */
#endif
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
    int64_t n, d;
    if (is_rational(e, &n, &d)) { *out = (double)n / (double)d; return true; }
    return false;
}

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool lg_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if `e` is exactly the symbol `name`. */
static bool lg_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

/* True if `e` is -Infinity, represented as Times[-1, Infinity]. */
static bool lg_is_neg_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!lg_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == -1 && lg_is_symbol(b, "Infinity"))
        return true;
    if (b->type == EXPR_INTEGER && b->data.integer == -1 && lg_is_symbol(a, "Infinity"))
        return true;
    return false;
}

/* True if `e` is a non-real directed infinity (e.g. I Infinity), represented
 * either as Times[Complex[..], Infinity] or DirectedInfinity[Complex[..]].
 * LogGamma maps all of these to ComplexInfinity. */
static bool lg_is_complex_infinity_dir(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr** a = e->data.function.args;
    size_t n = e->data.function.arg_count;
    if (lg_is_symbol(e->data.function.head, "Times") && n == 2) {
        bool has_inf = lg_is_symbol(a[0], "Infinity") || lg_is_symbol(a[1], "Infinity");
        bool has_cpx = (a[0]->type == EXPR_FUNCTION && lg_is_symbol(a[0]->data.function.head, "Complex")) ||
                       (a[1]->type == EXPR_FUNCTION && lg_is_symbol(a[1]->data.function.head, "Complex"));
        return has_inf && has_cpx;
    }
    if (lg_is_symbol(e->data.function.head, "DirectedInfinity") && n == 1) {
        return a[0]->type == EXPR_FUNCTION && lg_is_symbol(a[0]->data.function.head, "Complex");
    }
    return false;
}

/* ------------------------------------------------------------------ */
/* Exact path: integers and half-integers                             */
/* ------------------------------------------------------------------ */

/* Build Log[ Gamma[z] ] (evaluated), the exact value for positive integer and
 * positive half-integer z. Returns NULL if Gamma did not reduce. */
static Expr* lg_log_of_gamma(Expr* g_value) {
    Expr* arg = g_value;
    Expr* logc = expr_new_function(expr_new_symbol(SYM_Log), &arg, 1);
    return eval_and_free(logc);
}

/* Exact LogGamma for rational z with denominator 1 (integer) or 2
 * (half-integer). Returns NULL for anything else, leaving later paths to run. */
static Expr* lg_exact(Expr* arg) {
    int64_t n, d;
    if (!is_rational(arg, &n, &d)) return NULL;
    if (d != 1 && d != 2) return NULL;

    /* Non-positive integers are poles of Gamma; LogGamma diverges to Infinity. */
    if (d == 1 && n <= 0) return expr_new_symbol(SYM_Infinity);

    /* Evaluate the exact Gamma[z] (integer factorial or rational*Sqrt[Pi]). */
    Expr* zc = expr_copy(arg);
    Expr* g  = eval_and_free(expr_new_function(expr_new_symbol(SYM_Gamma), &zc, 1));
    if (!g) return NULL;
    if (g->type == EXPR_FUNCTION && lg_is_symbol(g->data.function.head, "Gamma")) {
        expr_free(g);   /* did not reduce -- stay symbolic */
        return NULL;
    }

    if (n > 0) {
        /* z > 0: Gamma(z) > 0, so LogGamma(z) = Log[Gamma(z)]. */
        return lg_log_of_gamma(g);
    }

    /* z < 0 half-integer (d == 2, n < 0): the analytic continuation differs
     * from Log[Gamma(z)] by a branch term. With Re = Log|Gamma(z)| and
     * Im = -Pi Ceiling[-z]:
     *   |Gamma(z)| = (-1)^Ceiling[-z] Gamma(z)  (sign alternates per interval)
     *   so posg = Gamma(z) when Ceiling[-z] even, else -Gamma(z) (positive). */
    int64_t p = -n;                       /* p > 0, z = -p/2 */
    int64_t c_ceil = (p + d - 1) / d;     /* Ceiling[p/d] = Ceiling[-z] */
    Expr* posg;
    if (c_ceil & 1) {
        /* Gamma(z) < 0: negate to obtain the positive magnitude. */
        posg = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), g }, 2));
    } else {
        posg = g;
    }
    Expr* re = lg_log_of_gamma(posg);
    if (!re) return NULL;

    /* im = Times[-c_ceil, I, Pi]. */
    Expr* im = expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ expr_new_integer(-c_ceil),
                             expr_new_symbol(SYM_I),
                             expr_new_symbol(SYM_Pi) }, 3);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                  (Expr*[]){ re, im }, 2));
}

/* ------------------------------------------------------------------ */
/* Machine complex path: Lanczos log-gamma                            */
/* ------------------------------------------------------------------ */

/* log(Gamma(z)) for complex z via the (g = 7, n = 9) Lanczos series. The log
 * of each Lanczos factor is summed directly, so this stays finite for large
 * |z| where Gamma itself overflows. Reflection handles Re(z) < 1/2. */
static double complex lg_lanczos(double complex z) {
    static const double g = 7.0;
    static const double c[9] = {
         0.99999999999980993,    676.5203681218851,    -1259.1392167224028,
       771.32342877765313,      -176.61502916214059,      12.507343278686905,
        -0.13857109526572012,      9.9843695780195716e-6,   1.5056327351493116e-7
    };
    if (creal(z) < 0.5) {
        /* Gamma(z) Gamma(1-z) = pi / sin(pi z); principal-branch log. */
        return clog(M_PI / csin(M_PI * z)) - lg_lanczos(1.0 - z);
    }
    double complex z1 = z - 1.0;
    double complex x = c[0];
    for (int i = 1; i < 9; i++) x += c[i] / (z1 + (double)i);
    double complex t = z1 + g + 0.5;
    return 0.5 * log(2.0 * M_PI) + (z1 + 0.5) * clog(t) - t + clog(x);
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ */
/* MPFR helpers (real)                                                */
/* ------------------------------------------------------------------ */

#define GRND MPFR_RNDN

/* ------------------------------------------------------------------ */
/* Complex-MPFR toolkit (pairs of mpfr_t; no MPC library available).  */
/* A `lcx` is one arbitrary-precision complex number; binary ops are  */
/* alias-safe (out may equal an input) and work at precision `p`.     */
/* This mirrors the `gcx` toolkit in gamma.c (kept local per the      */
/* established one-toolkit-per-module pattern, e.g. polygamma.c).     */
/* ------------------------------------------------------------------ */

typedef struct { mpfr_t re, im; } lcx;

static void lcx_init(lcx* z, mpfr_prec_t p) { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void lcx_clear(lcx* z)               { mpfr_clear(z->re);    mpfr_clear(z->im);    }
static void lcx_set(lcx* d, const lcx* s)   { mpfr_set(d->re, s->re, GRND); mpfr_set(d->im, s->im, GRND); }

static void lcx_add(lcx* out, const lcx* a, const lcx* b) {
    mpfr_add(out->re, a->re, b->re, GRND);
    mpfr_add(out->im, a->im, b->im, GRND);
}
static void lcx_sub(lcx* out, const lcx* a, const lcx* b) {
    mpfr_sub(out->re, a->re, b->re, GRND);
    mpfr_sub(out->im, a->im, b->im, GRND);
}

/* out = a * b. */
static void lcx_mul(lcx* out, const lcx* a, const lcx* b, mpfr_prec_t p) {
    mpfr_t ac, bd, ad, bc;
    mpfr_inits2(p, ac, bd, ad, bc, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, GRND);
    mpfr_mul(bd, a->im, b->im, GRND);
    mpfr_mul(ad, a->re, b->im, GRND);
    mpfr_mul(bc, a->im, b->re, GRND);
    mpfr_sub(out->re, ac, bd, GRND);
    mpfr_add(out->im, ad, bc, GRND);
    mpfr_clears(ac, bd, ad, bc, (mpfr_ptr)0);
}

/* out = a / b. */
static void lcx_div(lcx* out, const lcx* a, const lcx* b, mpfr_prec_t p) {
    mpfr_t den, t1, t2, nr, ni;
    mpfr_inits2(p, den, t1, t2, nr, ni, (mpfr_ptr)0);
    mpfr_mul(t1, b->re, b->re, GRND);
    mpfr_mul(t2, b->im, b->im, GRND);
    mpfr_add(den, t1, t2, GRND);                 /* |b|^2 */
    mpfr_mul(t1, a->re, b->re, GRND);
    mpfr_mul(t2, a->im, b->im, GRND);
    mpfr_add(nr, t1, t2, GRND);                  /* ac + bd */
    mpfr_mul(t1, a->im, b->re, GRND);
    mpfr_mul(t2, a->re, b->im, GRND);
    mpfr_sub(ni, t1, t2, GRND);                  /* bc - ad */
    mpfr_div(out->re, nr, den, GRND);
    mpfr_div(out->im, ni, den, GRND);
    mpfr_clears(den, t1, t2, nr, ni, (mpfr_ptr)0);
}

/* out = Log(z), principal branch: ln|z| + i Arg(z). */
static void lcx_log(lcx* out, const lcx* z, mpfr_prec_t p) {
    mpfr_t mag, ang;
    mpfr_inits2(p, mag, ang, (mpfr_ptr)0);
    mpfr_hypot(mag, z->re, z->im, GRND);
    mpfr_atan2(ang, z->im, z->re, GRND);
    mpfr_log(out->re, mag, GRND);
    mpfr_set(out->im, ang, GRND);
    mpfr_clears(mag, ang, (mpfr_ptr)0);
}

/* out = sin(z) = sin(a)cosh(b) + i cos(a)sinh(b). */
static void lcx_sin(lcx* out, const lcx* z, mpfr_prec_t p) {
    mpfr_t sa, ca, chb, shb;
    mpfr_inits2(p, sa, ca, chb, shb, (mpfr_ptr)0);
    mpfr_sin_cos(sa, ca, z->re, GRND);
    mpfr_sinh_cosh(shb, chb, z->im, GRND);
    mpfr_mul(out->re, sa, chb, GRND);
    mpfr_mul(out->im, ca, shb, GRND);
    mpfr_clears(sa, ca, chb, shb, (mpfr_ptr)0);
}

/* Fill be[1..K] with the even Bernoulli numbers B_2, B_4, ..., B_{2K} as exact
 * GMP rationals (the be[] entries must be mpq_init'd by the caller). Uses the
 * standard recurrence Sum_{j=0}^{m} C(m+1,j) B_j = 0 (m >= 1) over the full
 * B_0..B_{2K} sequence, keeping only the even-index values. */
static void lg_bernoulli_even(mpq_t* be, size_t K) {
    size_t Bn = 2 * K;
    mpq_t* B = (mpq_t*)malloc((Bn + 1) * sizeof(mpq_t));
    for (size_t i = 0; i <= Bn; i++) mpq_init(B[i]);
    mpq_set_ui(B[0], 1, 1);

    mpq_t s, t;
    mpq_inits(s, t, (mpq_ptr)0);
    mpz_t c;
    mpz_init(c);
    for (size_t m = 1; m <= Bn; m++) {
        mpq_set_ui(s, 0, 1);
        mpz_set_ui(c, 1);                         /* C(m+1, 0) = 1 */
        for (size_t j = 0; j < m; j++) {
            mpq_set_z(t, c);
            mpq_mul(t, t, B[j]);
            mpq_add(s, s, t);
            /* C(m+1, j+1) = C(m+1, j) * (m+1-j) / (j+1). */
            mpz_mul_ui(c, c, (unsigned long)(m + 1 - j));
            mpz_divexact_ui(c, c, (unsigned long)(j + 1));
        }
        mpq_set_ui(t, 1, (unsigned long)(m + 1)); /* 1/(m+1) */
        mpq_mul(B[m], s, t);
        mpq_neg(B[m], B[m]);                      /* B_m = -s/(m+1) */
    }
    mpz_clear(c);
    mpq_clears(s, t, (mpq_ptr)0);

    for (size_t k = 1; k <= K; k++) mpq_set(be[k], B[2 * k]);
    for (size_t i = 0; i <= Bn; i++) mpq_clear(B[i]);
    free(B);
}

/* LogGamma(z) for complex z into already-init'd `out` (precision wp), via the
 * Stirling asymptotic series with argument reduction. This is the analytic
 * (continuous-branch) log-gamma: its imaginary part grows without bound with
 * Im(z), which a single Log[Gamma[z]] (principal branch) cannot represent.
 *
 *   LogGamma(z) = LogGamma(w) - Sum_{j=0}^{M-1} Log(z+j),   w = z + M,
 *   LogGamma(w) = (w-1/2) Log(w) - w + (1/2) Log(2 pi)
 *                 + Sum_{k>=1} B_{2k} / (2k(2k-1) w^{2k-1}).
 *
 * z is shifted up by M so that Re(w) is large enough for the asymptotic series
 * to reach wp-bit accuracy before its terms start to diverge. Because Re of
 * every z+j is >= 1/2 > 0, each subtracted Log is a principal log with
 * argument in (-pi/2, pi/2), so the reduction is branch-continuous. Reflection
 * handles Re(z) < 1/2. */
static void lcx_loggamma(lcx* out, const lcx* z, mpfr_prec_t wp) {
    if (mpfr_cmp_d(z->re, 0.5) < 0) {
        /* LogGamma(z) = Log(pi/sin(pi z)) - LogGamma(1-z) (principal branch). */
        lcx omz, lg1, piz, spz, ratio, pic;
        lcx_init(&omz, wp); lcx_init(&lg1, wp); lcx_init(&piz, wp);
        lcx_init(&spz, wp); lcx_init(&ratio, wp); lcx_init(&pic, wp);
        mpfr_const_pi(pic.re, GRND); mpfr_set_ui(pic.im, 0, GRND);
        mpfr_ui_sub(omz.re, 1, z->re, GRND);     /* 1 - z */
        mpfr_neg(omz.im, z->im, GRND);
        lcx_loggamma(&lg1, &omz, wp);            /* LogGamma(1-z), Re >= 1/2 */
        lcx_mul(&piz, &pic, z, wp);              /* pi z */
        lcx_sin(&spz, &piz, wp);                 /* sin(pi z) */
        lcx_div(&ratio, &pic, &spz, wp);         /* pi / sin(pi z) */
        lcx_log(out, &ratio, wp);                /* Log(pi/sin(pi z)) */
        lcx_sub(out, out, &lg1);
        lcx_clear(&omz); lcx_clear(&lg1); lcx_clear(&piz);
        lcx_clear(&spz); lcx_clear(&ratio); lcx_clear(&pic);
        return;
    }

    /* K Bernoulli terms and a reduction target Wbound chosen so the optimally
     * truncated Stirling series (terms shrink while 2pi|w| > 2k) clears wp
     * bits: with Wbound ~ K/3 we have 2pi Wbound > 2K, and the minimal term
     * ~ exp(-2pi Wbound) << 2^-wp. */
    size_t K = (size_t)(0.5 * (double)wp) + 30;
    double Wbound = (double)K / 3.0 + 10.0;

    /* Number of integer shifts so Re(z) + M >= Wbound. */
    long M = 0;
    {
        mpfr_t rb;
        mpfr_init2(rb, wp);
        mpfr_set_d(rb, Wbound, GRND);
        if (mpfr_cmp(z->re, rb) < 0) {
            mpfr_sub(rb, rb, z->re, GRND);
            mpfr_ceil(rb, rb);
            M = (long)mpfr_get_si(rb, GRND);
            if (M < 0) M = 0;
        }
        mpfr_clear(rb);
    }

    /* sumlogs = Sum_{j=0}^{M-1} Log(z + j). */
    lcx sumlogs, zj, lzj;
    lcx_init(&sumlogs, wp); lcx_init(&zj, wp); lcx_init(&lzj, wp);
    mpfr_set_ui(sumlogs.re, 0, GRND); mpfr_set_ui(sumlogs.im, 0, GRND);
    for (long j = 0; j < M; j++) {
        mpfr_add_si(zj.re, z->re, j, GRND);
        mpfr_set(zj.im, z->im, GRND);
        lcx_log(&lzj, &zj, wp);
        lcx_add(&sumlogs, &sumlogs, &lzj);
    }

    /* w = z + M. */
    lcx w;
    lcx_init(&w, wp);
    mpfr_add_si(w.re, z->re, M, GRND);
    mpfr_set(w.im, z->im, GRND);

    /* out = (w-1/2) Log(w) - w + (1/2) Log(2 pi). */
    lcx half, logw;
    lcx_init(&half, wp); lcx_init(&logw, wp);
    mpfr_sub_d(half.re, w.re, 0.5, GRND);
    mpfr_set(half.im, w.im, GRND);
    lcx_log(&logw, &w, wp);
    lcx_mul(out, &half, &logw, wp);
    lcx_sub(out, out, &w);
    {
        mpfr_t twopi;
        mpfr_init2(twopi, wp);
        mpfr_const_pi(twopi, GRND);
        mpfr_mul_ui(twopi, twopi, 2, GRND);
        mpfr_log(twopi, twopi, GRND);
        mpfr_div_2ui(twopi, twopi, 1, GRND);     /* (1/2) Log(2 pi) */
        mpfr_add(out->re, out->re, twopi, GRND);
        mpfr_clear(twopi);
    }

    /* Stirling tail: Sum_k B_{2k}/(2k(2k-1)) * w^{-(2k-1)}.
     * p tracks w^{-(2k-1)} via p *= w^-2 each step. */
    mpq_t* be = (mpq_t*)malloc((K + 1) * sizeof(mpq_t));
    for (size_t i = 0; i <= K; i++) mpq_init(be[i]);
    lg_bernoulli_even(be, K);

    lcx invw, inv2, p, term;
    lcx_init(&invw, wp); lcx_init(&inv2, wp); lcx_init(&p, wp); lcx_init(&term, wp);
    {
        lcx one;
        lcx_init(&one, wp);
        mpfr_set_ui(one.re, 1, GRND); mpfr_set_ui(one.im, 0, GRND);
        lcx_div(&invw, &one, &w, wp);            /* 1/w */
        lcx_clear(&one);
    }
    lcx_mul(&inv2, &invw, &invw, wp);            /* 1/w^2 */
    lcx_set(&p, &invw);                          /* w^-1 (k=1) */

    mpfr_t coeff, tmag, prevmag, eps, scale;
    mpfr_inits2(wp, coeff, tmag, prevmag, eps, scale, (mpfr_ptr)0);
    mpfr_hypot(scale, out->re, out->im, GRND);   /* magnitude of leading part */
    if (mpfr_cmp_ui(scale, 1) < 0) mpfr_set_ui(scale, 1, GRND);
    mpfr_set_ui(eps, 1, GRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 8 ? wp - 8 : 1), GRND);
    mpfr_mul(eps, eps, scale, GRND);             /* relative tolerance */
    mpfr_set_inf(prevmag, 1);

    for (size_t k = 1; k <= K; k++) {
        unsigned long twok = (unsigned long)(2 * k);
        /* coeff = B_{2k} / (2k (2k-1)). */
        mpfr_set_q(coeff, be[k], GRND);
        mpfr_div_ui(coeff, coeff, twok, GRND);
        mpfr_div_ui(coeff, coeff, twok - 1, GRND);
        /* term = coeff * p. */
        mpfr_mul(term.re, p.re, coeff, GRND);
        mpfr_mul(term.im, p.im, coeff, GRND);
        mpfr_hypot(tmag, term.re, term.im, GRND);
        /* Asymptotic series: stop once terms start growing again. */
        if (mpfr_cmp(tmag, prevmag) > 0) break;
        lcx_add(out, out, &term);
        if (mpfr_cmp(tmag, eps) < 0) break;
        mpfr_set(prevmag, tmag, GRND);
        lcx_mul(&p, &p, &inv2, wp);              /* w^-(2k-1) -> w^-(2k+1) */
    }

    mpfr_clears(coeff, tmag, prevmag, eps, scale, (mpfr_ptr)0);
    for (size_t i = 0; i <= K; i++) mpq_clear(be[i]);
    free(be);

    /* out = LogGamma(w) - Sum_{j=0}^{M-1} Log(z+j). */
    lcx_sub(out, out, &sumlogs);

    lcx_clear(&invw); lcx_clear(&inv2); lcx_clear(&p); lcx_clear(&term);
    lcx_clear(&half); lcx_clear(&logw); lcx_clear(&w);
    lcx_clear(&sumlogs); lcx_clear(&zj); lcx_clear(&lzj);
}

/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool lg_set_mpfr(mpfr_t out, const Expr* e) {
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

/* Result precision (bits) for LogGamma[re + im I] under Mathematica
 * contagion: the minimum precision among the inexact parts, floored at
 * machine (53). See numeric_min_inexact_bits. */
static mpfr_prec_t lg_complex_prec(const Expr* re, const Expr* im) {
    long pr = numeric_min_inexact_bits(re);
    long pi = numeric_min_inexact_bits(im);
    long m  = (pr && pi) ? (pr < pi ? pr : pi) : (pr ? pr : pi);
    return m < 53 ? 53 : (mpfr_prec_t)m;
}

/* Build a Complex result from lcx components at the requested precision. */
static Expr* lg_complex_result(const mpfr_t re, const mpfr_t im, mpfr_prec_t out_prec) {
    Expr *rr, *ii;
    if (out_prec <= 53) {
        rr = expr_new_real(mpfr_get_d(re, GRND));
        ii = expr_new_real(mpfr_get_d(im, GRND));
    } else {
        rr = expr_new_mpfr_bits(out_prec);
        ii = expr_new_mpfr_bits(out_prec);
        mpfr_set(rr->data.mpfr, re, GRND);
        mpfr_set(ii->data.mpfr, im, GRND);
    }
    return make_complex(rr, ii);
}

/* Evaluate LogGamma(z) for a numeric complex `arg` at arbitrary precision.
 * Returns an EXPR_MPFR Complex/Real, ComplexInfinity at a pole, or NULL. */
static Expr* lg_mpfr_complex(Expr* re, Expr* im) {
    mpfr_prec_t out_prec = lg_complex_prec(re, im);
    mpfr_prec_t wp = out_prec + 64;              /* guard bits */

    lcx z, g;
    lcx_init(&z, wp); lcx_init(&g, wp);
    bool ok = lg_set_mpfr(z.re, re) && lg_set_mpfr(z.im, im);
    Expr* out = NULL;
    if (ok) {
        lcx_loggamma(&g, &z, wp);
        if (mpfr_inf_p(g.re) || mpfr_inf_p(g.im) ||
            mpfr_nan_p(g.re) || mpfr_nan_p(g.im)) {
            out = expr_new_symbol(SYM_ComplexInfinity);
        } else {
            out = lg_complex_result(g.re, g.im, out_prec);
        }
    }
    lcx_clear(&z); lcx_clear(&g);
    return out;
}

/* LogGamma for an arbitrary-precision real (EXPR_MPFR). */
static Expr* lg_mpfr_real(Expr* arg) {
    mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);

    /* Non-positive integer -> pole. */
    if (mpfr_integer_p(arg->data.mpfr) && mpfr_sgn(arg->data.mpfr) <= 0)
        return expr_new_symbol(SYM_ComplexInfinity);

    int sgn;
    mpfr_t r;
    mpfr_init2(r, prec);
    mpfr_lgamma(r, &sgn, arg->data.mpfr, GRND);  /* log|Gamma(z)| */
    if (mpfr_inf_p(r) || mpfr_nan_p(r)) {
        mpfr_clear(r);
        return expr_new_symbol(SYM_ComplexInfinity);
    }

    if (mpfr_sgn(arg->data.mpfr) > 0) {
        Expr* out = expr_new_mpfr_bits(prec);
        mpfr_set(out->data.mpfr, r, GRND);
        mpfr_clear(r);
        return out;
    }

    /* z < 0 (non-integer): complex result, Im = -Pi Ceiling[-z]. */
    Expr* reE = expr_new_mpfr_bits(prec);
    mpfr_set(reE->data.mpfr, r, GRND);
    mpfr_t nz, im;
    mpfr_init2(nz, prec); mpfr_init2(im, prec);
    mpfr_neg(nz, arg->data.mpfr, GRND);
    mpfr_ceil(nz, nz);                            /* Ceiling[-z] */
    mpfr_const_pi(im, GRND);
    mpfr_mul(im, im, nz, GRND);
    mpfr_neg(im, im, GRND);
    Expr* imE = expr_new_mpfr_bits(prec);
    mpfr_set(imE->data.mpfr, im, GRND);
    mpfr_clears(r, nz, im, (mpfr_ptr)0);
    return make_complex(reE, imE);
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* LogGamma[z]                                                        */
/* ------------------------------------------------------------------ */

static Expr* loggamma_one_arg(Expr* arg) {
    /* 1. Exact integer / half-integer (incl. the Infinity pole). */
    {
        Expr* exact = lg_exact(arg);
        if (exact) return exact;
    }

    /* 2. Symbolic infinities. */
    if (lg_is_symbol(arg, "Infinity"))        return expr_new_symbol(SYM_Infinity);
    if (lg_is_symbol(arg, "ComplexInfinity")) return expr_new_symbol(SYM_ComplexInfinity);
    if (lg_is_symbol(arg, "Indeterminate"))   return expr_new_symbol(SYM_Indeterminate);
    if (lg_is_neg_infinity(arg))              return expr_new_symbol(SYM_Indeterminate);
    if (lg_is_complex_infinity_dir(arg))      return expr_new_symbol(SYM_ComplexInfinity);

    /* 3. Machine real. */
    if (arg->type == EXPR_REAL) {
        double v = arg->data.real;
        if (v > 0.0) {
            double r = lgamma(v);
            if (isinf(r) || isnan(r)) return NULL; /* overflow: stay symbolic */
            return expr_new_real(r);
        }
        if (v == floor(v)) return expr_new_symbol(SYM_ComplexInfinity); /* pole */
        /* z < 0 non-integer: log|Gamma| + branch term -> complex. */
        double re = lgamma(v);
        double im = -M_PI * ceil(-v);
        return make_complex(expr_new_real(re), expr_new_real(im));
    }

#ifdef USE_MPFR
    /* 4. Arbitrary-precision real. */
    if (arg->type == EXPR_MPFR) return lg_mpfr_real(arg);
#endif

    /* 5. Complex argument (Complex[..] with an inexact part). */
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        bool inexact = lg_is_inexact(re) || lg_is_inexact(im);
#ifdef USE_MPFR
        if (inexact && (re->type == EXPR_MPFR || im->type == EXPR_MPFR)) {
            Expr* out = lg_mpfr_complex(re, im);
            if (out) return out;
        }
#endif
        double rr, ii;
        if (inexact && lg_to_double(re, &rr) && lg_to_double(im, &ii)) {
            double complex lg = lg_lanczos(rr + ii * I);
            if (cimag(lg) == 0.0) return expr_new_real(creal(lg));
            return make_complex(expr_new_real(creal(lg)), expr_new_real(cimag(lg)));
        }
    }

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

static Expr* loggamma_emit_argx(size_t argc) {
    fprintf(stderr,
            "LogGamma::argx: LogGamma called with %zu arguments; "
            "1 argument is expected.\n",
            argc);
    return NULL;
}

Expr* builtin_loggamma(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return loggamma_one_arg(args[0]);
    return loggamma_emit_argx(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void loggamma_init(void) {
    symtab_add_builtin("LogGamma", builtin_loggamma);
    symtab_get_def("LogGamma")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

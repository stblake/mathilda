/* Mathilda -- the Gamma function family.
 *
 *   Gamma[z]          Euler gamma function   Gamma(z) = Int_0^Inf t^(z-1) e^-t dt
 *   Gamma[a, z]       upper incomplete gamma Gamma(a,z) = Int_z^Inf t^(a-1) e^-t dt
 *   Gamma[a, z0, z1]  generalized incomplete = Gamma(a,z0) - Gamma(a,z1)
 *
 * Evaluation is layered so each kind of argument takes the cheapest exact
 * or fastest numeric route available:
 *
 *   exact integer / half-integer  ->  (z-1)! via the Factorial machinery
 *                                      (exact, BigInt, or rational*Sqrt[Pi])
 *   machine real        ->  libm   tgamma
 *   machine complex     ->  Lanczos approximation (double complex)
 *   arbitrary real      ->  MPFR   mpfr_gamma / mpfr_gamma_inc
 *   everything else     ->  stays symbolic (return NULL)
 *
 * Arbitrary-precision *complex* gamma is deliberately left symbolic: a
 * fixed-coefficient Lanczos series only carries ~15 correct digits, so
 * emitting it as an MPFR value would advertise a precision it does not
 * have. Reporting the input unevaluated is the honest behaviour.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "gamma.h"

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "arithmetic.h"   /* is_rational, make_rational, is_complex, make_complex */
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
/* Small numeric-coercion helpers                                     */
/* ------------------------------------------------------------------ */

/* Coerce an exact-or-real leaf to a double. Succeeds for Integer, Real,
 * BigInt and Rational; fails (returns false) for symbols, complex values,
 * MPFR values and anything else. */
static bool gamma_to_double(const Expr* e, double* out) {
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

/* True if `e` is an inexact numeric leaf (Real or MPFR): its presence is
 * what turns an otherwise-symbolic call into a numeric one. */
static bool gamma_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if `e` is exactly the symbol `name`. */
static bool gamma_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol, name) == 0;
}

/* True if `e` is -Infinity, represented as Times[-1, Infinity]. */
static bool gamma_is_neg_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!gamma_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == -1 && gamma_is_symbol(b, "Infinity"))
        return true;
    if (b->type == EXPR_INTEGER && b->data.integer == -1 && gamma_is_symbol(a, "Infinity"))
        return true;
    return false;
}

/* ------------------------------------------------------------------ */
/* Exact path: Gamma[z] = (z-1)! for integer / half-integer z         */
/* ------------------------------------------------------------------ */

/* Reuse the existing Factorial builtin (exact integers, BigInt, and the
 * Sqrt[Pi] half-integer rationals) by evaluating Factorial[z-1]. Only
 * called for is_rational(z) with denominator 1 or 2, the cases Factorial
 * actually closes. Returns NULL (and frees its scratch) if the Factorial
 * unexpectedly fails to reduce, leaving Gamma symbolic. */
static Expr* gamma_exact_via_factorial(int64_t n, int64_t d) {
    /* z = n/d, so z - 1 = (n - d)/d. */
    Expr* zm1 = make_rational(n - d, d);
    if (!zm1) return NULL;
    Expr* call = expr_new_function(expr_new_symbol("Factorial"), &zm1, 1);
    Expr* out = eval_and_free(call);
    /* If it came back still headed by Factorial, it did not reduce. */
    if (out && out->type == EXPR_FUNCTION &&
        gamma_is_symbol(out->data.function.head, "Factorial")) {
        expr_free(out);
        return NULL;
    }
    return out;
}

/* ------------------------------------------------------------------ */
/* Machine complex path: Lanczos approximation                        */
/* ------------------------------------------------------------------ */

/* Lanczos approximation (g = 7, n = 9), accurate to ~15 significant
 * digits across the complex plane. Uses the reflection formula for the
 * left half-plane where the series is ill-conditioned. */
static double complex gamma_lanczos(double complex z) {
    static const double g = 7.0;
    static const double c[9] = {
         0.99999999999980993,    676.5203681218851,    -1259.1392167224028,
       771.32342877765313,      -176.61502916214059,      12.507343278686905,
        -0.13857109526572012,      9.9843695780195716e-6,   1.5056327351493116e-7
    };
    if (creal(z) < 0.5) {
        /* Gamma(z) Gamma(1-z) = pi / sin(pi z). */
        return M_PI / (csin(M_PI * z) * gamma_lanczos(1.0 - z));
    }
    z -= 1.0;
    double complex x = c[0];
    for (int i = 1; i < 9; i++) x += c[i] / (z + (double)i);
    double complex t = z + g + 0.5;
    return csqrt(2.0 * M_PI) * cpow(t, z + 0.5) * cexp(-t) * x;
}

/* ------------------------------------------------------------------ */
/* MPFR helpers                                                       */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* Set an already-init2'd mpfr from an exact-or-real leaf. Succeeds for
 * Integer, Real, BigInt, MPFR and Rational; fails for complex / symbolic. */
static bool gamma_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, MPFR_RNDN); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          MPFR_RNDN); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        MPFR_RNDN); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          MPFR_RNDN); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, MPFR_RNDN);
        mpfr_div_si(out, out, (long)d, MPFR_RNDN);
        return true;
    }
    return false;
}

/* Working precision for a numeric Gamma: the largest precision among any
 * MPFR operands, else 53 bits (machine). */
static mpfr_prec_t gamma_work_prec(const Expr* a, const Expr* b) {
    mpfr_prec_t p = 53;
    if (a && a->type == EXPR_MPFR && mpfr_get_prec(a->data.mpfr) > p) p = mpfr_get_prec(a->data.mpfr);
    if (b && b->type == EXPR_MPFR && mpfr_get_prec(b->data.mpfr) > p) p = mpfr_get_prec(b->data.mpfr);
    return p;
}

/* ------------------------------------------------------------------ */
/* Complex-MPFR toolkit (pairs of mpfr_t; no MPC library available)   */
/*                                                                    */
/* A `gcx` is one arbitrary-precision complex number. All binary ops  */
/* are alias-safe (out may equal an input): inputs are read into      */
/* temporaries before any output component is written. Every op works */
/* at an explicit working precision `p`.                              */
/* ------------------------------------------------------------------ */

#define GRND MPFR_RNDN

typedef struct { mpfr_t re, im; } gcx;

static void gcx_init(gcx* z, mpfr_prec_t p)  { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void gcx_clear(gcx* z)                { mpfr_clear(z->re);    mpfr_clear(z->im);    }
static void gcx_set(gcx* d, const gcx* s)    { mpfr_set(d->re, s->re, GRND); mpfr_set(d->im, s->im, GRND); }

/* out = a + b, out = a - b. */
static void gcx_add(gcx* out, const gcx* a, const gcx* b) {
    mpfr_add(out->re, a->re, b->re, GRND);
    mpfr_add(out->im, a->im, b->im, GRND);
}
static void gcx_sub(gcx* out, const gcx* a, const gcx* b) {
    mpfr_sub(out->re, a->re, b->re, GRND);
    mpfr_sub(out->im, a->im, b->im, GRND);
}

/* |z| into mag (precision of mag). */
static void gcx_abs(mpfr_t mag, const gcx* z) { mpfr_hypot(mag, z->re, z->im, GRND); }

/* out = a * b. */
static void gcx_mul(gcx* out, const gcx* a, const gcx* b, mpfr_prec_t p) {
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
static void gcx_div(gcx* out, const gcx* a, const gcx* b, mpfr_prec_t p) {
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

/* out = exp(z). */
static void gcx_exp(gcx* out, const gcx* z, mpfr_prec_t p) {
    mpfr_t ea, c, s;
    mpfr_inits2(p, ea, c, s, (mpfr_ptr)0);
    mpfr_exp(ea, z->re, GRND);
    mpfr_sin_cos(s, c, z->im, GRND);
    mpfr_mul(out->re, ea, c, GRND);
    mpfr_mul(out->im, ea, s, GRND);
    mpfr_clears(ea, c, s, (mpfr_ptr)0);
}

/* out = Log(z), principal branch: (1/2)ln|z|^2 + i Arg(z). */
static void gcx_log(gcx* out, const gcx* z, mpfr_prec_t p) {
    mpfr_t mag, ang;
    mpfr_inits2(p, mag, ang, (mpfr_ptr)0);
    mpfr_hypot(mag, z->re, z->im, GRND);         /* |z|  */
    mpfr_atan2(ang, z->im, z->re, GRND);         /* Arg  */
    mpfr_log(out->re, mag, GRND);
    mpfr_set(out->im, ang, GRND);
    mpfr_clears(mag, ang, (mpfr_ptr)0);
}

/* out = sin(z) = sin(a)cosh(b) + i cos(a)sinh(b). */
static void gcx_sin(gcx* out, const gcx* z, mpfr_prec_t p) {
    mpfr_t sa, ca, chb, shb;
    mpfr_inits2(p, sa, ca, chb, shb, (mpfr_ptr)0);
    mpfr_sin_cos(sa, ca, z->re, GRND);
    mpfr_sinh_cosh(shb, chb, z->im, GRND);
    mpfr_mul(out->re, sa, chb, GRND);
    mpfr_mul(out->im, ca, shb, GRND);
    mpfr_clears(sa, ca, chb, shb, (mpfr_ptr)0);
}

/* out = base^expo (principal branch) = exp(expo * Log(base)). */
static void gcx_pow(gcx* out, const gcx* base, const gcx* expo, mpfr_prec_t p) {
    gcx lg, prod;
    gcx_init(&lg, p); gcx_init(&prod, p);
    gcx_log(&lg, base, p);
    gcx_mul(&prod, expo, &lg, p);
    gcx_exp(out, &prod, p);
    gcx_clear(&lg); gcx_clear(&prod);
}

/* ------------------------------------------------------------------ */
/* Arbitrary-precision complex Gamma via Spouge's approximation       */
/*                                                                    */
/* Unlike a fixed-coefficient Lanczos series, Spouge's coefficients   */
/* are computable at runtime for any target accuracy, so this honours */
/* the requested precision. Reflection handles Re(z) < 1/2.           */
/*                                                                    */
/*   Gamma(z) = (z-1+N)^(z-1/2) e^-(z-1+N) [ c_0 + Sum_k c_k/(z-1+k) ] */
/*   c_0 = sqrt(2 pi),                                                 */
/*   c_k = (-1)^(k-1) (N-k)^(k-1/2) e^(N-k) / (k-1)!,  k = 1..N-1.     */
/* ------------------------------------------------------------------ */

/* Number of Spouge terms for working precision wp bits. The truncation
 * error is ~ (2 pi)^-N, so N ln(2 pi) >= wp ln 2; 0.45 wp clears that
 * with margin. */
static size_t gamma_spouge_terms(mpfr_prec_t wp) {
    size_t n = (size_t)(0.45 * (double)wp) + 4;
    return n;
}

/* Gamma(z) for complex z into already-init'd `out` (precision wp). */
static void gcx_gamma(gcx* out, const gcx* z, mpfr_prec_t wp) {
    /* Reflection for the left half-plane, where Spouge is ill-conditioned. */
    if (mpfr_cmp_d(z->re, 0.5) < 0) {
        gcx omz, g1, piz, spz, denom, pic;
        gcx_init(&omz, wp); gcx_init(&g1, wp); gcx_init(&piz, wp);
        gcx_init(&spz, wp); gcx_init(&denom, wp); gcx_init(&pic, wp);
        mpfr_const_pi(pic.re, GRND); mpfr_set_ui(pic.im, 0, GRND);
        mpfr_ui_sub(omz.re, 1, z->re, GRND);     /* 1 - z */
        mpfr_neg(omz.im, z->im, GRND);
        gcx_gamma(&g1, &omz, wp);                /* Gamma(1-z), Re >= 1/2 */
        gcx_mul(&piz, &pic, z, wp);              /* pi z */
        gcx_sin(&spz, &piz, wp);                 /* sin(pi z) */
        gcx_mul(&denom, &spz, &g1, wp);          /* sin(pi z) Gamma(1-z) */
        gcx_div(out, &pic, &denom, wp);          /* pi / denom */
        gcx_clear(&omz); gcx_clear(&g1); gcx_clear(&piz);
        gcx_clear(&spz); gcx_clear(&denom); gcx_clear(&pic);
        return;
    }

    size_t N = gamma_spouge_terms(wp);

    /* Accumulate S = c_0 + Sum_{k=1}^{N-1} c_k / (z-1+k). */
    gcx S, term, denom, num;
    gcx_init(&S, wp); gcx_init(&term, wp); gcx_init(&denom, wp); gcx_init(&num, wp);

    mpfr_t twopi;
    mpfr_init2(twopi, wp);
    mpfr_const_pi(twopi, GRND);
    mpfr_mul_ui(twopi, twopi, 2, GRND);
    mpfr_sqrt(S.re, twopi, GRND);                /* c_0 = sqrt(2 pi) */
    mpfr_set_ui(S.im, 0, GRND);
    mpfr_clear(twopi);

    mpfr_t ck, base, expo, ex, fact;
    mpfr_inits2(wp, ck, base, expo, ex, fact, (mpfr_ptr)0);
    mpfr_set_ui(fact, 1, GRND);                  /* (k-1)! starting at k=1: 0! = 1 */

    for (size_t k = 1; k < N; k++) {
        unsigned long ak = (unsigned long)(N - k); /* N - k > 0 */
        /* c_k = (-1)^(k-1) ak^(k-1/2) e^ak / (k-1)! */
        mpfr_set_ui(base, ak, GRND);
        mpfr_set_d(expo, (double)k - 0.5, GRND);
        mpfr_pow(ck, base, expo, GRND);          /* ak^(k-1/2) */
        mpfr_exp(ex, base, GRND);                /* e^ak */
        mpfr_mul(ck, ck, ex, GRND);
        mpfr_div(ck, ck, fact, GRND);            /* / (k-1)! */
        if ((k - 1) & 1) mpfr_neg(ck, ck, GRND); /* (-1)^(k-1) */

        /* term = (ck, 0) / (z - 1 + k) */
        mpfr_set(denom.re, z->re, GRND);
        mpfr_add_si(denom.re, denom.re, (long)k - 1, GRND);
        mpfr_set(denom.im, z->im, GRND);
        mpfr_set(num.re, ck, GRND);
        mpfr_set_ui(num.im, 0, GRND);
        gcx_div(&term, &num, &denom, wp);
        gcx_add(&S, &S, &term);

        mpfr_mul_ui(fact, fact, (unsigned long)k, GRND); /* (k-1)! -> k! */
    }
    mpfr_clears(ck, base, expo, ex, fact, (mpfr_ptr)0);

    /* pref = (z-1+N)^(z-1/2) * e^-(z-1+N) */
    gcx bb, expo_c, p1, p2;
    gcx_init(&bb, wp); gcx_init(&expo_c, wp); gcx_init(&p1, wp); gcx_init(&p2, wp);
    mpfr_add_si(bb.re, z->re, (long)N - 1, GRND);
    mpfr_set(bb.im, z->im, GRND);
    mpfr_sub_d(expo_c.re, z->re, 0.5, GRND);
    mpfr_set(expo_c.im, z->im, GRND);
    gcx_pow(&p1, &bb, &expo_c, wp);              /* (z-1+N)^(z-1/2) */
    mpfr_neg(bb.re, bb.re, GRND);                /* -(z-1+N) */
    mpfr_neg(bb.im, bb.im, GRND);
    gcx_exp(&p2, &bb, wp);                        /* e^-(z-1+N) */
    gcx_mul(out, &p1, &p2, wp);
    gcx_mul(out, out, &S, wp);

    gcx_clear(&bb); gcx_clear(&expo_c); gcx_clear(&p1); gcx_clear(&p2);
    gcx_clear(&S); gcx_clear(&term); gcx_clear(&denom); gcx_clear(&num);
}

/* Largest MPFR precision among the parts of a complex value (re, im),
 * else 53 (machine). */
static mpfr_prec_t gamma_complex_prec(const Expr* re, const Expr* im) {
    return gamma_work_prec(re, im);
}

/* Build a Complex result from gcx components: machine inputs (out_prec <= 53)
 * yield Real parts, arbitrary precision yields MPFR parts -- mirroring the
 * real Gamma/incomplete paths. make_complex drops a zero imaginary part. */
static Expr* gamma_complex_result(const mpfr_t re, const mpfr_t im, mpfr_prec_t out_prec) {
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

/* Evaluate Gamma(z) for a numeric complex `arg` (Complex[..] with at least
 * one MPFR part) at arbitrary precision. Returns an EXPR_MPFR Complex, or
 * ComplexInfinity at a pole, or NULL if the parts are not numeric. */
static Expr* gamma_mpfr_complex(Expr* re, Expr* im) {
    mpfr_prec_t out_prec = gamma_complex_prec(re, im);
    mpfr_prec_t wp = out_prec + 64;              /* guard bits */

    gcx z, g;
    gcx_init(&z, wp); gcx_init(&g, wp);
    bool ok = gamma_set_mpfr(z.re, re) && gamma_set_mpfr(z.im, im);
    Expr* out = NULL;
    if (ok) {
        gcx_gamma(&g, &z, wp);
        if (mpfr_inf_p(g.re) || mpfr_inf_p(g.im) ||
            mpfr_nan_p(g.re) || mpfr_nan_p(g.im)) {
            out = expr_new_symbol("ComplexInfinity");
        } else {
            out = gamma_complex_result(g.re, g.im, out_prec);
        }
    }
    gcx_clear(&z); gcx_clear(&g);
    return out;
}

/* ------------------------------------------------------------------ */
/* Arbitrary-precision complex incomplete gamma  Gamma(a, z)          */
/*                                                                    */
/* Two convergent representations, selected by region:                */
/*   - lower-series  gamma(a,z) = z^a e^-z Sum_n z^n/((a)_{n+1}),      */
/*     then Gamma(a,z) = Gamma(a) - gamma(a,z)     [Re(z) < Re(a)+1]   */
/*   - Lentz continued fraction for Gamma(a,z) directly  [otherwise].  */
/* Both finish with the e^-z z^a prefactor (principal branch).        */
/* ------------------------------------------------------------------ */

/* eps = 2^-(wp-4), the relative tolerance for the iterations. */
static void gamma_set_eps(mpfr_t eps, mpfr_prec_t wp) {
    mpfr_set_ui(eps, 1, GRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 8 ? wp - 4 : 1), GRND);
}

/* out *= e^-z z^a (the common incomplete-gamma prefactor). */
static void gcx_mul_incprefactor(gcx* out, const gcx* a, const gcx* z, mpfr_prec_t wp) {
    gcx za, ez, nz;
    gcx_init(&za, wp); gcx_init(&ez, wp); gcx_init(&nz, wp);
    gcx_pow(&za, z, a, wp);                       /* z^a */
    mpfr_neg(nz.re, z->re, GRND); mpfr_neg(nz.im, z->im, GRND);
    gcx_exp(&ez, &nz, wp);                        /* e^-z */
    gcx_mul(out, out, &za, wp);
    gcx_mul(out, out, &ez, wp);
    gcx_clear(&za); gcx_clear(&ez); gcx_clear(&nz);
}

/* Lower-incomplete-gamma series. Returns false if it does not converge. */
static bool gcx_lower_series(gcx* out, const gcx* a, const gcx* z, mpfr_prec_t wp) {
    gcx S, t, ap, tmp, one;
    gcx_init(&S, wp); gcx_init(&t, wp); gcx_init(&ap, wp);
    gcx_init(&tmp, wp); gcx_init(&one, wp);
    mpfr_set_ui(one.re, 1, GRND); mpfr_set_ui(one.im, 0, GRND);
    gcx_div(&t, &one, a, wp);                     /* t_0 = 1/a */
    gcx_set(&S, &t);
    mpfr_t eps, mag, smag;
    mpfr_inits2(wp, eps, mag, smag, (mpfr_ptr)0);
    gamma_set_eps(eps, wp);
    bool ok = false;
    for (unsigned long n = 1; n < 200000; n++) {
        mpfr_add_ui(ap.re, a->re, n, GRND);       /* a + n */
        mpfr_set(ap.im, a->im, GRND);
        gcx_mul(&tmp, &t, z, wp);
        gcx_div(&t, &tmp, &ap, wp);               /* t_n = t_{n-1} z/(a+n) */
        gcx_add(&S, &S, &t);
        gcx_abs(mag, &t); gcx_abs(smag, &S);
        mpfr_mul(smag, smag, eps, GRND);
        if (mpfr_cmp(mag, smag) < 0) { ok = true; break; }
    }
    if (ok) { gcx_set(out, &S); gcx_mul_incprefactor(out, a, z, wp); }
    gcx_clear(&S); gcx_clear(&t); gcx_clear(&ap);
    gcx_clear(&tmp); gcx_clear(&one);
    mpfr_clears(eps, mag, smag, (mpfr_ptr)0);
    return ok;
}

/* Upper-incomplete-gamma via modified Lentz continued fraction.
 * Valid for Re(z) > 0. Returns false if it does not converge. */
static bool gcx_upper_cf(gcx* out, const gcx* a, const gcx* z, mpfr_prec_t wp) {
    gcx b, c, d, an, del, tmp, ia, tinyc;
    gcx_init(&b, wp); gcx_init(&c, wp); gcx_init(&d, wp); gcx_init(&an, wp);
    gcx_init(&del, wp); gcx_init(&tmp, wp); gcx_init(&ia, wp); gcx_init(&tinyc, wp);
    mpfr_t tiny, eps, mag, t;
    mpfr_inits2(wp, tiny, eps, mag, t, (mpfr_ptr)0);
    mpfr_set_ui(tiny, 1, GRND);
    mpfr_div_2ui(tiny, tiny, (unsigned long)wp + 16, GRND);  /* avoids 0 pivots */
    gamma_set_eps(eps, wp);
    mpfr_set(tinyc.re, tiny, GRND); mpfr_set_ui(tinyc.im, 0, GRND);

    /* b = z + 1 - a;  c = 1/tiny;  d = 1/b;  h(=out) = d */
    mpfr_add_ui(b.re, z->re, 1, GRND); mpfr_sub(b.re, b.re, a->re, GRND);
    mpfr_sub(b.im, z->im, a->im, GRND);
    mpfr_ui_div(c.re, 1, tiny, GRND); mpfr_set_ui(c.im, 0, GRND);
    {
        gcx one; gcx_init(&one, wp);
        mpfr_set_ui(one.re, 1, GRND); mpfr_set_ui(one.im, 0, GRND);
        gcx_div(&d, &one, &b, wp);
        gcx_clear(&one);
    }
    gcx_set(out, &d);

    bool ok = false;
    for (unsigned long i = 1; i < 200000; i++) {
        /* an = -i (i - a) */
        mpfr_ui_sub(ia.re, i, a->re, GRND);       /* i - a */
        mpfr_neg(ia.im, a->im, GRND);
        mpfr_mul_ui(an.re, ia.re, i, GRND);       /* i (i-a) */
        mpfr_mul_ui(an.im, ia.im, i, GRND);
        mpfr_neg(an.re, an.re, GRND);             /* negate */
        mpfr_neg(an.im, an.im, GRND);
        /* b += 2 */
        mpfr_add_ui(b.re, b.re, 2, GRND);
        /* d = an*d + b */
        gcx_mul(&tmp, &an, &d, wp);
        gcx_add(&d, &tmp, &b);
        gcx_abs(mag, &d);
        if (mpfr_cmp(mag, tiny) < 0) gcx_set(&d, &tinyc);
        /* c = b + an/c */
        gcx_div(&tmp, &an, &c, wp);
        gcx_add(&c, &b, &tmp);
        gcx_abs(mag, &c);
        if (mpfr_cmp(mag, tiny) < 0) gcx_set(&c, &tinyc);
        /* d = 1/d;  del = d*c;  h *= del */
        {
            gcx one; gcx_init(&one, wp);
            mpfr_set_ui(one.re, 1, GRND); mpfr_set_ui(one.im, 0, GRND);
            gcx_div(&d, &one, &d, wp);
            gcx_clear(&one);
        }
        gcx_mul(&del, &d, &c, wp);
        gcx_mul(out, out, &del, wp);
        /* |del - 1| < eps ? */
        mpfr_sub_ui(t, del.re, 1, GRND);
        mpfr_hypot(mag, t, del.im, GRND);
        if (mpfr_cmp(mag, eps) < 0) { ok = true; break; }
    }
    if (ok) gcx_mul_incprefactor(out, a, z, wp);
    gcx_clear(&b); gcx_clear(&c); gcx_clear(&d); gcx_clear(&an);
    gcx_clear(&del); gcx_clear(&tmp); gcx_clear(&ia); gcx_clear(&tinyc);
    mpfr_clears(tiny, eps, mag, t, (mpfr_ptr)0);
    return ok;
}

/* Gamma(a, z) for numeric complex a/z. Returns false on non-convergence
 * (caller then leaves the expression symbolic). */
static bool gcx_inc_gamma(gcx* out, const gcx* a, const gcx* z, mpfr_prec_t wp) {
    /* Region split: series where Re(z) < Re(a)+1 (and Re(z) small enough that
     * the series is well-behaved), continued fraction otherwise. */
    mpfr_t thresh;
    mpfr_init2(thresh, wp);
    mpfr_add_ui(thresh, a->re, 1, GRND);
    bool use_series = (mpfr_cmp(z->re, thresh) < 0);
    mpfr_clear(thresh);

    if (use_series) {
        /* Gamma(a,z) = Gamma(a) - gamma(a,z). */
        gcx lower, ga;
        gcx_init(&lower, wp); gcx_init(&ga, wp);
        bool ok = gcx_lower_series(&lower, a, z, wp);
        if (ok) {
            gcx_gamma(&ga, a, wp);
            if (mpfr_inf_p(ga.re) || mpfr_inf_p(ga.im) ||
                mpfr_nan_p(ga.re) || mpfr_nan_p(ga.im)) {
                ok = false;                       /* Gamma(a) pole: series route invalid */
            } else {
                gcx_sub(out, &ga, &lower);
            }
        }
        gcx_clear(&lower); gcx_clear(&ga);
        if (ok) return true;
        /* fall through to CF if the series route did not work */
    }
    return gcx_upper_cf(out, a, z, wp);
}

/* Evaluate Gamma(a, z) where a and/or z are numeric complex values (with at
 * least one inexact part). Returns an EXPR_MPFR Complex / Real, or NULL. */
static Expr* gamma_mpfr_inc_complex(Expr* ar, Expr* ai, Expr* zr, Expr* zi,
                                    mpfr_prec_t out_prec) {
    mpfr_prec_t wp = out_prec + 64;
    gcx a, z, g;
    gcx_init(&a, wp); gcx_init(&z, wp); gcx_init(&g, wp);
    bool ok = gamma_set_mpfr(a.re, ar) && gamma_set_mpfr(a.im, ai) &&
              gamma_set_mpfr(z.re, zr) && gamma_set_mpfr(z.im, zi);
    Expr* out = NULL;
    if (ok && gcx_inc_gamma(&g, &a, &z, wp) &&
        !mpfr_nan_p(g.re) && !mpfr_nan_p(g.im) &&
        !mpfr_inf_p(g.re) && !mpfr_inf_p(g.im)) {
        out = gamma_complex_result(g.re, g.im, out_prec);
    }
    gcx_clear(&a); gcx_clear(&z); gcx_clear(&g);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Gamma[z]  (the Euler gamma function)                               */
/* ------------------------------------------------------------------ */

static Expr* gamma_one_arg(Expr* arg) {
    /* 1. Exact integer / half-integer -> (z-1)!. */
    int64_t n, d;
    if (is_rational(arg, &n, &d) && (d == 1 || d == 2)) {
        Expr* exact = gamma_exact_via_factorial(n, d);
        if (exact) return exact;
    }

    /* 2. Symbolic infinities. */
    if (gamma_is_symbol(arg, "Infinity"))        return expr_new_symbol("Infinity");
    if (gamma_is_symbol(arg, "ComplexInfinity")) return expr_new_symbol("ComplexInfinity");
    if (gamma_is_symbol(arg, "Indeterminate"))   return expr_new_symbol("Indeterminate");
    if (gamma_is_neg_infinity(arg))              return expr_new_symbol("Indeterminate");

    /* 3. Machine real. */
    if (arg->type == EXPR_REAL) {
        double v = arg->data.real;
        double r = tgamma(v);
        if (isnan(r)) return expr_new_symbol("ComplexInfinity"); /* pole at <=0 integer */
        if (isinf(r)) {
            if (v <= 0.0) return expr_new_symbol("ComplexInfinity");
            return NULL; /* overflow for large positive z: stay symbolic */
        }
        return expr_new_real(r);
    }

#ifdef USE_MPFR
    /* 4. Arbitrary-precision real. */
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        Expr* out = expr_new_mpfr_bits(prec);
        mpfr_gamma(out->data.mpfr, arg->data.mpfr, MPFR_RNDN);
        if (mpfr_inf_p(out->data.mpfr)) { expr_free(out); return expr_new_symbol("ComplexInfinity"); }
        if (mpfr_nan_p(out->data.mpfr)) { expr_free(out); return expr_new_symbol("ComplexInfinity"); }
        return out;
    }
#endif

    /* 5. Complex argument (Complex[..] with an inexact part). */
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        bool inexact = gamma_is_inexact(re) || gamma_is_inexact(im);
#ifdef USE_MPFR
        /* Arbitrary-precision complex: route through Spouge so the result
         * carries the requested precision. */
        if (inexact && (re->type == EXPR_MPFR || im->type == EXPR_MPFR)) {
            Expr* out = gamma_mpfr_complex(re, im);
            if (out) return out;
        }
#endif
        /* Machine complex: fast double Lanczos. */
        double rr, ii;
        if (inexact && gamma_to_double(re, &rr) && gamma_to_double(im, &ii)) {
            double complex g = gamma_lanczos(rr + ii * I);
            if (cimag(g) == 0.0) return expr_new_real(creal(g));
            return make_complex(expr_new_real(creal(g)), expr_new_real(cimag(g)));
        }
    }

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Exact path: Gamma[n, z] for positive integer n (finite closed form)*/
/* ------------------------------------------------------------------ */

/* Largest positive integer first argument for which we expand the
 * incomplete gamma into its (degree n-1) polynomial form. Beyond this
 * the expansion is more noise than help, so we leave it symbolic. */
#define GAMMA_INT_EXPAND_CAP 1000

/* For positive integer n, Gamma[n, z] = (n-1)! e^{-z} Sum_{k=0}^{n-1} z^k/k!
 *                                      = e^{-z} Sum_{k=0}^{n-1} c_k z^k,
 * with the integer coefficient c_k = (n-1)!/k!. Builds the symbolic
 * expression and evaluates it (so exact z stays exact: Gamma[2,3] -> 4/E^3,
 * symbolic z reduces: Gamma[2,x] -> (1+x)/E^x). Returns NULL if n is out of
 * range. The coefficients are computed with GMP so large n stay exact. */
static Expr* gamma_incomplete_int(int64_t n, Expr* z) {
    if (n < 1 || n > GAMMA_INT_EXPAND_CAP) return NULL;

    /* Integer coefficients c_k = (n-1)!/k!, built from the top down:
     *   c_{n-1} = 1,  c_{k} = c_{k+1} * (k+1). */
    size_t terms = (size_t)n;
    Expr** sum = (Expr**)malloc(terms * sizeof(Expr*));
    if (!sum) return NULL;

    mpz_t c;
    mpz_init_set_ui(c, 1);                       /* c_{n-1} = 1 */
    for (int64_t k = n - 1; k >= 0; k--) {
        Expr* coeff = expr_bigint_normalize(expr_new_bigint_from_mpz(c));
        Expr* term;
        if (k == 0) {
            term = coeff;                        /* z^0 = 1 */
        } else {
            Expr* zk = expr_new_function(expr_new_symbol("Power"),
                          (Expr*[]){ expr_copy(z), expr_new_integer(k) }, 2);
            term = expr_new_function(expr_new_symbol("Times"),
                          (Expr*[]){ coeff, zk }, 2);
        }
        sum[k] = term;
        if (k > 0) mpz_mul_ui(c, c, (unsigned long)k); /* c_{k-1} = c_k * k */
    }
    mpz_clear(c);

    Expr* poly = (terms == 1)
        ? sum[0]
        : expr_new_function(expr_new_symbol("Plus"), sum, terms);
    free(sum);

    /* e^{-z} = Exp[-z]. */
    Expr* nz  = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(z) }, 2));
    Expr* exp = eval_and_free(expr_new_function(expr_new_symbol("Exp"), &nz, 1));

    return eval_and_free(expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ exp, poly }, 2));
}

#ifdef USE_MPFR
/* Decompose a numeric value into borrowed (re, im) parts. A plain real value
 * yields (value, shared zero); a Complex[..] yields its two parts. */
static void gamma_decompose(Expr* v, Expr* zero, Expr** re, Expr** im) {
    Expr *r, *i;
    if (is_complex(v, &r, &i)) { *re = r; *im = i; }
    else { *re = v; *im = zero; }
}
#endif

/* ------------------------------------------------------------------ */
/* Gamma[a, z]  (upper incomplete gamma)                              */
/* ------------------------------------------------------------------ */

static Expr* gamma_two_arg(Expr* a, Expr* z) {
    /* 1. Exact rewrites that hold for any z / a. */
    if (z->type == EXPR_INTEGER && z->data.integer == 0) {
        /* Gamma[a, 0] = Gamma[a]. */
        Expr* ga = expr_copy(a);
        return eval_and_free(expr_new_function(expr_new_symbol("Gamma"), &ga, 1));
    }
    if (gamma_is_symbol(z, "Infinity")) return expr_new_integer(0); /* Gamma[a, Inf] = 0 */

    /* Positive integer a: finite closed form (covers Gamma[1, z] = E^-z).
     * Only for symbolic / exact z -- inexact z falls through to the fast
     * numeric path below rather than building a float polynomial. */
    if (a->type == EXPR_INTEGER && a->data.integer >= 1 && !gamma_is_inexact(z)) {
        Expr* closed = gamma_incomplete_int(a->data.integer, z);
        if (closed) return closed;
    }

#ifdef USE_MPFR
    /* 2. Numeric real path -- needs at least one inexact operand so we do
     *    not silently turn exact Gamma[2, 3] into a float. */
    if (gamma_is_inexact(a) || gamma_is_inexact(z)) {
        mpfr_prec_t prec = gamma_work_prec(a, z);
        mpfr_t av, zv, rv;
        mpfr_init2(av, prec); mpfr_init2(zv, prec); mpfr_init2(rv, prec);
        bool ok = gamma_set_mpfr(av, a) && gamma_set_mpfr(zv, z);
        Expr* out = NULL;
        if (ok) {
            mpfr_gamma_inc(rv, av, zv, MPFR_RNDN);
            if (mpfr_nan_p(rv)) {
                out = NULL;
            } else if (a->type == EXPR_MPFR || z->type == EXPR_MPFR) {
                out = expr_new_mpfr_copy(rv);            /* arbitrary precision */
            } else {
                out = expr_new_real(mpfr_get_d(rv, MPFR_RNDN)); /* machine */
            }
        }
        mpfr_clear(av); mpfr_clear(zv); mpfr_clear(rv);
        if (out) return out;
    }

    /* 3. Numeric complex path: a and/or z is complex-valued and at least one
     *    part is inexact (so exact complex like Gamma[3/2, I] stays symbolic). */
    {
        Expr* zero = expr_new_integer(0);
        Expr *are, *aim, *zre, *zim;
        gamma_decompose(a, zero, &are, &aim);
        gamma_decompose(z, zero, &zre, &zim);
        bool any_complex  = (aim != zero) || (zim != zero);
        bool any_inexact  = gamma_is_inexact(are) || gamma_is_inexact(aim) ||
                            gamma_is_inexact(zre) || gamma_is_inexact(zim);
        if (any_complex && any_inexact) {
            mpfr_prec_t p1 = gamma_work_prec(are, aim);
            mpfr_prec_t p2 = gamma_work_prec(zre, zim);
            mpfr_prec_t p  = p1 > p2 ? p1 : p2;
            Expr* out = gamma_mpfr_inc_complex(are, aim, zre, zim, p);
            if (out) { expr_free(zero); return out; }
        }
        expr_free(zero);
    }
#else
    (void)gamma_is_inexact;
#endif

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Gamma[a, z0, z1] = Gamma[a, z0] - Gamma[a, z1]                      */
/* ------------------------------------------------------------------ */

static Expr* gamma_three_arg(Expr* a, Expr* z0, Expr* z1) {
    Expr* g0 = expr_new_function(expr_new_symbol("Gamma"),
                                 (Expr*[]){ expr_copy(a), expr_copy(z0) }, 2);
    Expr* g1 = expr_new_function(expr_new_symbol("Gamma"),
                                 (Expr*[]){ expr_copy(a), expr_copy(z1) }, 2);
    Expr* diff = expr_new_function(expr_new_symbol("Subtract"),
                                   (Expr*[]){ g0, g1 }, 2);
    return eval_and_free(diff);
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

Expr* builtin_gamma(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return gamma_one_arg(args[0]);
    if (argc == 2) return gamma_two_arg(args[0], args[1]);
    if (argc == 3) return gamma_three_arg(args[0], args[1], args[2]);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void gamma_init(void) {
    symtab_add_builtin("Gamma", builtin_gamma);
    symtab_get_def("Gamma")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

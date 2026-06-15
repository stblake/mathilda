/* Mathilda -- the Airy function Ai.
 *
 *   AiryAi[z]   Airy function Ai(z), the solution of  y'' = z y  that tends to
 *               zero as z -> +Infinity along the real axis. Ai is an *entire*
 *               function of z (no branch cuts).
 *
 * Evaluation is layered so each kind of argument takes the most accurate and
 * cheapest route:
 *
 *   exact special values   ->  AiryAi[0] = 1/(3^(2/3) Gamma[2/3]),
 *                              AiryAi[+-Infinity] = 0
 *   machine real           ->  unified complex-MPFR core at 53 bits, real part
 *   arbitrary real         ->  unified complex-MPFR core at mpfr_get_prec bits
 *   complex (any precision) ->  unified complex-MPFR core, Complex[..] result
 *   everything else        ->  stays symbolic (return NULL)
 *
 * The unified core `airy_ai_core` evaluates Ai(z) and Ai'(z) together in a
 * file-local complex-MPFR toolkit (`acx`, pairs of mpfr_t -- no MPC library is
 * available; this mirrors the `ecx`/`pcx`/`gcx` toolkits in erf.c/polylog.c/
 * gamma.c). It routes between two algorithms on r = |z| and the requested
 * output precision P:
 *
 *   - Maclaurin series (small/moderate |z|). From Ai'' = z Ai the Taylor
 *     coefficients satisfy a_0 = Ai(0), a_1 = Ai'(0), a_2 = 0 and
 *     a_n = a_{n-3} / (n (n-1)) for n >= 3. The partial sums reach magnitude
 *     ~exp((2/3) r^{3/2}) before cancelling for complex / negative arguments,
 *     so the core adds  (2/3) r^{3/2} / ln2  guard bits to absorb that exactly.
 *
 *   - Asymptotic series (large |z|), DLMF 9.7.5/9.7.6. With zeta = (2/3) z^{3/2}
 *         Ai(z)  ~ exp(-zeta)/(2 sqrt(pi) z^{1/4}) Sum (-1)^k u_k / zeta^k,
 *         Ai'(z) ~ -z^{1/4} exp(-zeta)/(2 sqrt(pi)) Sum (-1)^k v_k / zeta^k,
 *     summed to the optimal (smallest-term) truncation. The single series is
 *     accurate for |arg z| <= 2 pi / 3; closer to the negative real axis
 *     (2 pi / 3 < |arg z| <= pi) the core uses the connection relation
 *     (DLMF 9.2.12)
 *         Ai(z) = -[ w Ai(w z) + conj(w) Ai(conj(w) z)],  w = e^{2 pi i / 3},
 *     which maps the argument into two points with |arg| <= 2 pi / 3 where the
 *     direct series is accurate; the oscillation on the negative real axis then
 *     emerges naturally from the sum of the two rotated evaluations.
 *
 * D[AiryAi[z], z] = AiryAiPrime[z] (see calculus/deriv.c); the Maclaurin series
 * at 0 is produced by the generic Taylor-via-D path once AiryAi[0] / AiryAiPrime[0]
 * have closed-form values.
 *
 * AiryAiPrime[z] = Ai'(z) is a full numeric evaluator in its own right: because
 * `airy_ai_core` returns Ai(z) and Ai'(z) together, AiryAiPrime reuses the very
 * same Maclaurin / asymptotic / connection machinery and simply selects the
 * derivative component. Its exact values are AiryAiPrime[0] = -1/(3^(1/3) Gamma[1/3])
 * and AiryAiPrime[+Infinity] = 0; at -Infinity Ai' has no limit (oscillation with
 * growing ~|z|^(1/4) amplitude) and is left unevaluated.
 *
 * Attributes (both heads): Listable, NumericFunction, Protected, ReadProtected.
 */
#include "airyai.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arithmetic.h"   /* is_rational, make_complex, is_complex */
#include "numeric.h"      /* numeric_min_inexact_bits */
#include "attr.h"
#include "eval.h"         /* eval_and_free */
#include "expr.h"
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* M_PI / M_LN2 are POSIX, not C99 -- provide fallbacks (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

/* True if `e` is exactly the symbol `name`. */
static bool ai_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol, name) == 0;
}

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool ai_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if `e` is -Infinity, represented as Times[-1, Infinity]. */
static bool ai_is_neg_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!ai_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == -1 && ai_is_symbol(b, "Infinity"))
        return true;
    if (b->type == EXPR_INTEGER && b->data.integer == -1 && ai_is_symbol(a, "Infinity"))
        return true;
    return false;
}

/* Build the exact value AiryAi[0] = 1 / (3^(2/3) Gamma[2/3]) and evaluate it. */
static Expr* ai_value_at_zero(void) {
    /* 3^(2/3) */
    Expr* p = expr_new_function(expr_new_symbol("Power"),
                  (Expr*[]){ expr_new_integer(3), make_rational(2, 3) }, 2);
    /* Gamma[2/3] */
    Expr* g = expr_new_function(expr_new_symbol("Gamma"),
                  (Expr*[]){ make_rational(2, 3) }, 1);
    /* 3^(2/3) Gamma[2/3] */
    Expr* denom = expr_new_function(expr_new_symbol("Times"),
                  (Expr*[]){ p, g }, 2);
    /* 1 / (...) = Power[denom, -1] */
    Expr* inv = expr_new_function(expr_new_symbol("Power"),
                  (Expr*[]){ denom, expr_new_integer(-1) }, 2);
    return eval_and_free(inv);
}

/* Build the exact value AiryAiPrime[0] = -1 / (3^(1/3) Gamma[1/3]). */
static Expr* ai_prime_value_at_zero(void) {
    Expr* p = expr_new_function(expr_new_symbol("Power"),
                  (Expr*[]){ expr_new_integer(3), make_rational(1, 3) }, 2);
    Expr* g = expr_new_function(expr_new_symbol("Gamma"),
                  (Expr*[]){ make_rational(1, 3) }, 1);
    Expr* denom = expr_new_function(expr_new_symbol("Times"),
                  (Expr*[]){ p, g }, 2);
    Expr* inv = expr_new_function(expr_new_symbol("Power"),
                  (Expr*[]){ denom, expr_new_integer(-1) }, 2);
    Expr* neg = expr_new_function(expr_new_symbol("Times"),
                  (Expr*[]){ expr_new_integer(-1), inv }, 2);
    return eval_and_free(neg);
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ */
/* Complex-MPFR toolkit (pairs of mpfr_t; no MPC library available).  */
/* Each op runs at an explicit working precision and is alias-safe.    */
/* ------------------------------------------------------------------ */

#define ARND MPFR_RNDN

typedef struct { mpfr_t re, im; } acx;

static void acx_init(acx* z, mpfr_prec_t p) { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void acx_clear(acx* z)               { mpfr_clear(z->re);    mpfr_clear(z->im);    }
static void acx_set(acx* d, const acx* s)   { mpfr_set(d->re, s->re, ARND); mpfr_set(d->im, s->im, ARND); }

static void acx_add(acx* out, const acx* a, const acx* b) {
    mpfr_add(out->re, a->re, b->re, ARND);
    mpfr_add(out->im, a->im, b->im, ARND);
}

/* out = a * b. */
static void acx_mul(acx* out, const acx* a, const acx* b, mpfr_prec_t p) {
    mpfr_t ac, bd, ad, bc;
    mpfr_inits2(p, ac, bd, ad, bc, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, ARND);
    mpfr_mul(bd, a->im, b->im, ARND);
    mpfr_mul(ad, a->re, b->im, ARND);
    mpfr_mul(bc, a->im, b->re, ARND);
    mpfr_sub(out->re, ac, bd, ARND);
    mpfr_add(out->im, ad, bc, ARND);
    mpfr_clears(ac, bd, ad, bc, (mpfr_ptr)0);
}

/* out = a / b. */
static void acx_div(acx* out, const acx* a, const acx* b, mpfr_prec_t p) {
    mpfr_t ac, bd, ad, bc, den;
    mpfr_inits2(p, ac, bd, ad, bc, den, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, ARND);
    mpfr_mul(bd, a->im, b->im, ARND);
    mpfr_mul(ad, a->im, b->re, ARND);   /* note: a->im b->re */
    mpfr_mul(bc, a->re, b->im, ARND);   /*       a->re b->im */
    mpfr_mul(den, b->re, b->re, ARND);
    mpfr_fma(den, b->im, b->im, den, ARND);     /* |b|^2 */
    mpfr_add(ac, ac, bd, ARND);                 /* re num = ac + bd */
    mpfr_sub(ad, ad, bc, ARND);                 /* im num = ad - bc */
    mpfr_div(out->re, ac, den, ARND);
    mpfr_div(out->im, ad, den, ARND);
    mpfr_clears(ac, bd, ad, bc, den, (mpfr_ptr)0);
}

/* |z| into mag. */
static void acx_abs(mpfr_t mag, const acx* z) { mpfr_hypot(mag, z->re, z->im, ARND); }

/* arg(z) (principal, in (-pi, pi]) into out. */
static void acx_arg(mpfr_t out, const acx* z) { mpfr_atan2(out, z->im, z->re, ARND); }

/* out = exp(z) = e^re (cos im + i sin im). */
static void acx_exp(acx* out, const acx* z, mpfr_prec_t p) {
    mpfr_t ea, c, s;
    mpfr_inits2(p, ea, c, s, (mpfr_ptr)0);
    mpfr_exp(ea, z->re, ARND);
    mpfr_sin_cos(s, c, z->im, ARND);
    mpfr_mul(out->re, ea, c, ARND);
    mpfr_mul(out->im, ea, s, ARND);
    mpfr_clears(ea, c, s, (mpfr_ptr)0);
}

/* out = z^e (principal branch) for a real exponent e:
 *   out = |z|^e (cos(e arg z) + i sin(e arg z)). */
static void acx_pow_d(acx* out, const acx* z, double e, mpfr_prec_t p) {
    mpfr_t r, th, c, s, ex;
    mpfr_inits2(p, r, th, c, s, ex, (mpfr_ptr)0);
    acx_abs(r, z);
    acx_arg(th, z);
    mpfr_set_d(ex, e, ARND);
    mpfr_pow(r, r, ex, ARND);       /* |z|^e */
    mpfr_mul_d(th, th, e, ARND);    /* e arg z */
    mpfr_sin_cos(s, c, th, ARND);
    mpfr_mul(out->re, r, c, ARND);
    mpfr_mul(out->im, r, s, ARND);
    mpfr_clears(r, th, c, s, ex, (mpfr_ptr)0);
}

/* out = z scaled by the real mpfr s (alias-safe). */
static void acx_scale(acx* out, const acx* z, const mpfr_t s) {
    mpfr_mul(out->re, z->re, s, ARND);
    mpfr_mul(out->im, z->im, s, ARND);
}

/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool ai_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, ARND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          ARND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        ARND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          ARND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, ARND);
        mpfr_div_si(out, out, (long)d, ARND);
        return true;
    }
    return false;
}

/* The two Maclaurin constants, computed at precision p:
 *   a0 = Ai(0)  =  1 / (3^(2/3) Gamma[2/3]),
 *   a1 = Ai'(0) = -1 / (3^(1/3) Gamma[1/3]). */
static void airy_consts(mpfr_t a0, mpfr_t a1, mpfr_prec_t p) {
    mpfr_t three, e, g, frac;
    mpfr_inits2(p, three, e, g, frac, (mpfr_ptr)0);
    mpfr_set_ui(three, 3, ARND);

    /* a0 = 1 / (3^(2/3) Gamma(2/3)).  Form 2/3 *exactly* in MPFR -- passing
     * 2.0/3.0 as a double would cap the constant at ~53 bits. */
    mpfr_set_ui(frac, 2, ARND); mpfr_div_ui(frac, frac, 3, ARND);   /* 2/3 */
    mpfr_pow(e, three, frac, ARND);       /* 3^(2/3) */
    mpfr_gamma(g, frac, ARND);            /* Gamma(2/3) */
    mpfr_mul(a0, e, g, ARND);
    mpfr_ui_div(a0, 1, a0, ARND);

    /* a1 = -1 / (3^(1/3) Gamma(1/3)). */
    mpfr_set_ui(frac, 1, ARND); mpfr_div_ui(frac, frac, 3, ARND);   /* 1/3 */
    mpfr_pow(e, three, frac, ARND);       /* 3^(1/3) */
    mpfr_gamma(g, frac, ARND);            /* Gamma(1/3) */
    mpfr_mul(a1, e, g, ARND);
    mpfr_ui_div(a1, 1, a1, ARND);
    mpfr_neg(a1, a1, ARND);

    mpfr_clears(three, e, g, frac, (mpfr_ptr)0);
}

/* ------------------------------------------------------------------ */
/* Maclaurin series (small/moderate |z|), accurate everywhere.        */
/*   Ai(z)  = Sum_n a_n z^n,  Ai'(z) = Sum_n n a_n z^{n-1},           */
/*   a_0 = Ai(0), a_1 = Ai'(0), a_2 = 0, a_n = a_{n-3}/(n(n-1)).      */
/* Outputs Ai, Aip are init'd by the caller at precision wp.          */
/* ------------------------------------------------------------------ */
static void airy_maclaurin(acx* Ai, acx* Aip, const acx* z, mpfr_prec_t wp) {
    mpfr_t a[3];        /* sliding window: a[n%3] holds the current coefficient */
    mpfr_init2(a[0], wp); mpfr_init2(a[1], wp); mpfr_init2(a[2], wp);
    airy_consts(a[0], a[1], wp);
    mpfr_set_ui(a[2], 0, ARND);

    acx zpow, zprev, term;          /* z^n, z^{n-1}, a_n z^n */
    acx_init(&zpow, wp); acx_init(&zprev, wp); acx_init(&term, wp);
    mpfr_set_ui(zpow.re, 1, ARND); mpfr_set_ui(zpow.im, 0, ARND);   /* z^0 */
    mpfr_set_ui(zprev.re, 0, ARND); mpfr_set_ui(zprev.im, 0, ARND);

    mpfr_set_ui(Ai->re, 0, ARND);  mpfr_set_ui(Ai->im, 0, ARND);
    mpfr_set_ui(Aip->re, 0, ARND); mpfr_set_ui(Aip->im, 0, ARND);

    mpfr_t r, mag, smag, eps, coef;
    mpfr_inits2(wp, r, mag, smag, eps, coef, (mpfr_ptr)0);
    acx_abs(r, z);
    double rd = mpfr_get_d(r, ARND);
    mpfr_set_ui(eps, 1, ARND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 8 ? wp - 4 : 1), ARND);  /* 2^-(wp-4) */

    double peak = 2.0 * pow(rd, 1.5) + 8.0;
    unsigned long cap = (unsigned long)(4.0 * pow(rd, 1.5)) + (unsigned long)wp + 1000;

    for (unsigned long n = 0; n <= cap; n++) {
        if (n >= 3) {
            /* a_n = a_{n-3} / (n (n-1)); a[n%3] currently holds a_{n-3}. */
            mpfr_div_ui(a[n % 3], a[n % 3], n, ARND);
            mpfr_div_ui(a[n % 3], a[n % 3], n - 1, ARND);
        }
        mpfr_set(coef, a[n % 3], ARND);                 /* a_n */

        /* term = a_n * z^n; Ai += term. */
        acx_scale(&term, &zpow, coef);
        acx_add(Ai, Ai, &term);

        /* Ai' += (n a_n) * z^{n-1}. */
        if (n >= 1) {
            mpfr_mul_ui(smag, coef, n, ARND);           /* n a_n (smag reused as scratch) */
            acx_scale(&term, &zprev, smag);
            acx_add(Aip, Aip, &term);
        }

        /* Convergence: past the term peak and the Ai term is negligible.
         * Skip the structural zeros (every n == 2 mod 3 has a_n = 0) -- they
         * would spuriously satisfy the test and truncate the series early. */
        if ((double)n > peak && mpfr_sgn(coef) != 0) {
            acx_scale(&term, &zpow, coef);
            acx_abs(mag, &term);
            acx_abs(smag, Ai);
            mpfr_mul(smag, smag, eps, ARND);
            if (mpfr_cmp(mag, smag) < 0) break;
        }

        /* Advance powers: zprev = z^n, zpow = z^{n+1}. */
        acx_set(&zprev, &zpow);
        acx_mul(&zpow, &zpow, z, wp);
    }

    mpfr_clears(r, mag, smag, eps, coef, (mpfr_ptr)0);
    mpfr_clear(a[0]); mpfr_clear(a[1]); mpfr_clear(a[2]);
    acx_clear(&zpow); acx_clear(&zprev); acx_clear(&term);
}

/* ------------------------------------------------------------------ */
/* Asymptotic series (DLMF 9.7.5/9.7.6), valid for |arg z| <= 2pi/3.  */
/*   zeta = (2/3) z^{3/2}                                              */
/*   Ai(z)  ~ exp(-zeta)/(2 sqrt(pi) z^{1/4}) Sum (-1)^k u_k/zeta^k    */
/*   Ai'(z) ~ -z^{1/4} exp(-zeta)/(2 sqrt(pi)) Sum (-1)^k v_k/zeta^k   */
/* Summed to the optimal (smallest-term) truncation.                  */
/* ------------------------------------------------------------------ */
static void airy_asymp_direct(acx* Ai, acx* Aip, const acx* z, mpfr_prec_t wp) {
    acx zeta, nzeta, z14, ez, invz, invpow, term, sumA, sumAp;
    acx_init(&zeta, wp); acx_init(&nzeta, wp); acx_init(&z14, wp);
    acx_init(&ez, wp); acx_init(&invz, wp); acx_init(&invpow, wp);
    acx_init(&term, wp); acx_init(&sumA, wp); acx_init(&sumAp, wp);

    /* zeta = (2/3) z^{3/2}, z14 = z^{1/4}.  Scale by 2/3 exactly in MPFR. */
    acx_pow_d(&zeta, z, 1.5, wp);
    mpfr_mul_ui(zeta.re, zeta.re, 2, ARND); mpfr_div_ui(zeta.re, zeta.re, 3, ARND);
    mpfr_mul_ui(zeta.im, zeta.im, 2, ARND); mpfr_div_ui(zeta.im, zeta.im, 3, ARND);
    acx_pow_d(&z14, z, 0.25, wp);

    /* exp(-zeta). */
    mpfr_neg(nzeta.re, zeta.re, ARND);
    mpfr_neg(nzeta.im, zeta.im, ARND);
    acx_exp(&ez, &nzeta, wp);

    /* invz = 1/zeta. */
    {
        acx one; acx_init(&one, wp);
        mpfr_set_ui(one.re, 1, ARND); mpfr_set_ui(one.im, 0, ARND);
        acx_div(&invz, &one, &zeta, wp);
        acx_clear(&one);
    }

    /* k = 0 terms. */
    mpfr_set_ui(sumA.re, 1, ARND);  mpfr_set_ui(sumA.im, 0, ARND);
    mpfr_set_ui(sumAp.re, 1, ARND); mpfr_set_ui(sumAp.im, 0, ARND);
    mpfr_set_ui(invpow.re, 1, ARND); mpfr_set_ui(invpow.im, 0, ARND);   /* zeta^0 */

    mpfr_t u, v, coef, mag, prevmag;
    mpfr_inits2(wp, u, v, coef, mag, prevmag, (mpfr_ptr)0);
    mpfr_set_ui(u, 1, ARND);
    mpfr_set_inf(prevmag, 1);

    for (unsigned long k = 1; k < 100000; k++) {
        /* u_k = (6k-5)(6k-3)(6k-1)/((2k-1) 216 k) u_{k-1}  (DLMF 9.7.2). */
        mpfr_mul_ui(u, u, 6 * k - 5, ARND);
        mpfr_mul_ui(u, u, 6 * k - 3, ARND);
        mpfr_mul_ui(u, u, 6 * k - 1, ARND);
        mpfr_div_ui(u, u, 216, ARND);
        mpfr_div_ui(u, u, k, ARND);
        mpfr_div_ui(u, u, 2 * k - 1, ARND);
        /* v_k = -(6k+1)/(6k-1) u_k. */
        mpfr_mul_ui(v, u, 6 * k + 1, ARND);
        mpfr_div_ui(v, v, 6 * k - 1, ARND);
        mpfr_neg(v, v, ARND);

        /* invpow = zeta^{-k}. */
        acx_mul(&invpow, &invpow, &invz, wp);

        /* termA = (-1)^k u_k zeta^{-k}; magnitude drives optimal truncation. */
        mpfr_set(coef, u, ARND);
        if (k & 1) mpfr_neg(coef, coef, ARND);
        acx_scale(&term, &invpow, coef);
        acx_abs(mag, &term);
        if (mpfr_cmp(mag, prevmag) >= 0) break;     /* smallest term reached */
        acx_add(&sumA, &sumA, &term);

        /* termAp = (-1)^k v_k zeta^{-k}. */
        mpfr_set(coef, v, ARND);
        if (k & 1) mpfr_neg(coef, coef, ARND);
        acx_scale(&term, &invpow, coef);
        acx_add(&sumAp, &sumAp, &term);

        mpfr_set(prevmag, mag, ARND);
    }

    /* 2 sqrt(pi). */
    mpfr_t twosqrtpi;
    mpfr_init2(twosqrtpi, wp);
    mpfr_const_pi(twosqrtpi, ARND);
    mpfr_sqrt(twosqrtpi, twosqrtpi, ARND);
    mpfr_mul_ui(twosqrtpi, twosqrtpi, 2, ARND);

    /* Ai = exp(-zeta)/(2 sqrt(pi) z^{1/4}) * sumA. */
    {
        acx denom, prefA;
        acx_init(&denom, wp); acx_init(&prefA, wp);
        acx_scale(&denom, &z14, twosqrtpi);     /* 2 sqrt(pi) z^{1/4} */
        acx_div(&prefA, &ez, &denom, wp);
        acx_mul(Ai, &prefA, &sumA, wp);
        acx_clear(&denom); acx_clear(&prefA);
    }
    /* Ai' = -z^{1/4} exp(-zeta)/(2 sqrt(pi)) * sumAp. */
    {
        acx prefAp;
        acx_init(&prefAp, wp);
        acx_mul(&prefAp, &z14, &ez, wp);
        mpfr_div(prefAp.re, prefAp.re, twosqrtpi, ARND);
        mpfr_div(prefAp.im, prefAp.im, twosqrtpi, ARND);
        mpfr_neg(prefAp.re, prefAp.re, ARND);
        mpfr_neg(prefAp.im, prefAp.im, ARND);
        acx_mul(Aip, &prefAp, &sumAp, wp);
        acx_clear(&prefAp);
    }

    mpfr_clear(twosqrtpi);
    mpfr_clears(u, v, coef, mag, prevmag, (mpfr_ptr)0);
    acx_clear(&zeta); acx_clear(&nzeta); acx_clear(&z14);
    acx_clear(&ez); acx_clear(&invz); acx_clear(&invpow);
    acx_clear(&term); acx_clear(&sumA); acx_clear(&sumAp);
}

/* ------------------------------------------------------------------ */
/* Unified core: Ai(z) and Ai'(z) at output precision P bits.         */
/* Picks Maclaurin (small |z|) or asymptotic (large |z|), and uses    */
/* the DLMF 9.2.12 connection relation near the negative real axis.   */
/* Outputs Ai, Aip are init'd by the caller (any precision).          */
/* ------------------------------------------------------------------ */
static void airy_ai_core(const acx* z, acx* Ai, acx* Aip, mpfr_prec_t P) {
    mpfr_t rm; mpfr_init2(rm, (P < 64 ? 64 : P));
    acx_abs(rm, z);
    double rd = mpfr_get_d(rm, ARND);
    mpfr_clear(rm);

    /* Smallest |z| at which the asymptotic series can reach P bits:
     * optimal-truncation error ~ exp(-(4/3) r^{3/2}); a 1.3 margin keeps
     * us clear of the seam. */
    double rmin = pow(0.75 * (double)P * M_LN2, 2.0 / 3.0);
    if (rmin < 1.0) rmin = 1.0;

    if (rd > 1.3 * rmin) {
        /* Asymptotic region. */
        mpfr_t am; mpfr_init2(am, (P < 64 ? 64 : P));
        acx_arg(am, z);
        double arg = mpfr_get_d(am, ARND);
        mpfr_clear(am);

        const double TWO_PI_3 = 2.0 * M_PI / 3.0;
        if (fabs(arg) <= TWO_PI_3 + 1e-12) {
            mpfr_prec_t wp = P + 64;
            acx ai, aip; acx_init(&ai, wp); acx_init(&aip, wp);
            airy_asymp_direct(&ai, &aip, z, wp);
            mpfr_set(Ai->re, ai.re, ARND);  mpfr_set(Ai->im, ai.im, ARND);
            mpfr_set(Aip->re, aip.re, ARND); mpfr_set(Aip->im, aip.im, ARND);
            acx_clear(&ai); acx_clear(&aip);
        } else {
            /* Connection relation Ai(z) = -[w Ai(wz) + wbar Ai(wbar z)],
             *   Ai'(z) = -[wbar Ai'(wz) + w Ai'(wbar z)],  w = e^{2pi i/3}.
             * Both wz and wbar z have |arg| <= 2pi/3 (direct series valid),
             * and the recessive rotated term carries no large cancellation. */
            mpfr_prec_t wp = P + 96;
            acx w, wb, z1, z2, ai1, aip1, ai2, aip2, t1, t2;
            acx_init(&w, wp); acx_init(&wb, wp);
            acx_init(&z1, wp); acx_init(&z2, wp);
            acx_init(&ai1, wp); acx_init(&aip1, wp);
            acx_init(&ai2, wp); acx_init(&aip2, wp);
            acx_init(&t1, wp); acx_init(&t2, wp);

            /* w = -1/2 + i sqrt(3)/2; wbar = conj(w). */
            mpfr_set_d(w.re, -0.5, ARND);
            mpfr_sqrt_ui(w.im, 3, ARND); mpfr_div_2ui(w.im, w.im, 1, ARND);
            mpfr_set(wb.re, w.re, ARND); mpfr_neg(wb.im, w.im, ARND);

            acx_mul(&z1, &w, z, wp);
            acx_mul(&z2, &wb, z, wp);
            airy_asymp_direct(&ai1, &aip1, &z1, wp);
            airy_asymp_direct(&ai2, &aip2, &z2, wp);

            /* Ai = -(w ai1 + wbar ai2). */
            acx_mul(&t1, &w, &ai1, wp);
            acx_mul(&t2, &wb, &ai2, wp);
            acx_add(&t1, &t1, &t2);
            mpfr_neg(Ai->re, t1.re, ARND); mpfr_neg(Ai->im, t1.im, ARND);

            /* Ai' = -(wbar aip1 + w aip2). */
            acx_mul(&t1, &wb, &aip1, wp);
            acx_mul(&t2, &w, &aip2, wp);
            acx_add(&t1, &t1, &t2);
            mpfr_neg(Aip->re, t1.re, ARND); mpfr_neg(Aip->im, t1.im, ARND);

            acx_clear(&w); acx_clear(&wb); acx_clear(&z1); acx_clear(&z2);
            acx_clear(&ai1); acx_clear(&aip1); acx_clear(&ai2); acx_clear(&aip2);
            acx_clear(&t1); acx_clear(&t2);
        }
    } else {
        /* Maclaurin region. Guard absorbs the ~exp((2/3) r^{3/2}) partial-sum
         * cancellation; capped so a pathological argument cannot blow up. */
        long guard = 64 + (long)((2.0 / 3.0) * pow(rd, 1.5) / M_LN2);
        mpfr_prec_t wp = P + (mpfr_prec_t)guard;
        if (wp > P + 300000) wp = P + 300000;
        acx ai, aip; acx_init(&ai, wp); acx_init(&aip, wp);
        airy_maclaurin(&ai, &aip, z, wp);
        mpfr_set(Ai->re, ai.re, ARND);  mpfr_set(Ai->im, ai.im, ARND);
        mpfr_set(Aip->re, aip.re, ARND); mpfr_set(Aip->im, aip.im, ARND);
        acx_clear(&ai); acx_clear(&aip);
    }
}

/* Build a numeric Expr from a real mpfr value at out_prec bits: a machine
 * double when out_prec <= 53 (promoting to MPFR only if the magnitude
 * overflows the double range), otherwise an MPFR leaf. */
static Expr* ai_real_result(const mpfr_t v, mpfr_prec_t out_prec) {
    if (out_prec <= 53) {
        double d = mpfr_get_d(v, ARND);
        if (isinf(d) && !mpfr_inf_p(v))   /* overflowed double: keep precision */
            return expr_new_mpfr_copy(v);
        return expr_new_real(d);
    }
    Expr* out = expr_new_mpfr_bits(out_prec);
    mpfr_set(out->data.mpfr, v, ARND);
    return out;
}

/* Build a complex result Complex[re, im] at out_prec bits (Real parts for
 * machine precision, MPFR parts above). If a machine-precision result
 * overflows the double range (|Ai| can reach ~10^374 for large imaginary z),
 * both parts are promoted to MPFR so the value is retained, matching
 * Mathematica (where such a result is no longer a machine number). */
static Expr* ai_complex_result(const mpfr_t re, const mpfr_t im, mpfr_prec_t out_prec) {
    Expr *rr, *ii;
    if (out_prec <= 53) {
        double dr = mpfr_get_d(re, ARND), di = mpfr_get_d(im, ARND);
        bool overflow = (isinf(dr) && !mpfr_inf_p(re)) || (isinf(di) && !mpfr_inf_p(im));
        if (overflow) {
            rr = expr_new_mpfr_copy(re);
            ii = expr_new_mpfr_copy(im);
            return make_complex(rr, ii);
        }
        rr = expr_new_real(dr);
        ii = expr_new_real(di);
    } else {
        rr = expr_new_mpfr_bits(out_prec);
        ii = expr_new_mpfr_bits(out_prec);
        mpfr_set(rr->data.mpfr, re, ARND);
        mpfr_set(ii->data.mpfr, im, ARND);
    }
    return make_complex(rr, ii);
}

/* AiryAi (prime=false) or AiryAiPrime (prime=true) for a real numeric leaf
 * (Real / MPFR) at out_prec bits. The core computes Ai and Ai' together; we
 * select the requested component for the result. */
static Expr* ai_eval_real(const Expr* arg, mpfr_prec_t out_prec, bool prime) {
    mpfr_prec_t wp = (out_prec < 64 ? 64 : out_prec);
    acx z, Ai, Aip;
    acx_init(&z, wp); acx_init(&Ai, wp); acx_init(&Aip, wp);
    Expr* out = NULL;
    if (ai_set_mpfr(z.re, arg)) {
        mpfr_set_ui(z.im, 0, ARND);
        airy_ai_core(&z, &Ai, &Aip, out_prec);
        const acx* R = prime ? &Aip : &Ai;
        if (!mpfr_nan_p(R->re) && !mpfr_inf_p(R->re))
            out = ai_real_result(R->re, out_prec);
    }
    acx_clear(&z); acx_clear(&Ai); acx_clear(&Aip);
    return out;
}

/* AiryAi (prime=false) or AiryAiPrime (prime=true) for a numeric
 * Complex[re, im] at out_prec bits. */
static Expr* ai_eval_complex(const Expr* re, const Expr* im, mpfr_prec_t out_prec,
                             bool prime) {
    mpfr_prec_t wp = (out_prec < 64 ? 64 : out_prec);
    acx z, Ai, Aip;
    acx_init(&z, wp); acx_init(&Ai, wp); acx_init(&Aip, wp);
    Expr* out = NULL;
    if (ai_set_mpfr(z.re, re) && ai_set_mpfr(z.im, im)) {
        airy_ai_core(&z, &Ai, &Aip, out_prec);
        const acx* R = prime ? &Aip : &Ai;
        if (!mpfr_nan_p(R->re) && !mpfr_nan_p(R->im) &&
            !mpfr_inf_p(R->re) && !mpfr_inf_p(R->im))
            out = ai_complex_result(R->re, R->im, out_prec);
    }
    acx_clear(&z); acx_clear(&Ai); acx_clear(&Aip);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* AiryAi[z]                                                          */
/* ------------------------------------------------------------------ */

static Expr* airyai_one_arg(Expr* arg) {
    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return ai_value_at_zero();                 /* AiryAi[0] = 1/(3^(2/3) Gamma[2/3]) */
    if (ai_is_symbol(arg, "Infinity"))   return expr_new_integer(0);  /* Ai(+inf) = 0 */
    if (ai_is_neg_infinity(arg))          return expr_new_integer(0);  /* Ai(-inf) = 0 */
    if (ai_is_symbol(arg, "Indeterminate")) return expr_new_symbol("Indeterminate");

#ifdef USE_MPFR
    /* 2. Machine real. */
    if (arg->type == EXPR_REAL)
        return ai_eval_real(arg, 53, false);

    /* 3. Arbitrary-precision real. */
    if (arg->type == EXPR_MPFR)
        return ai_eval_real(arg, mpfr_get_prec(arg->data.mpfr), false);

    /* 4. Complex argument (Complex[..] with an inexact part). */
    {
        Expr *re, *im;
        if (is_complex(arg, &re, &im) && (ai_is_inexact(re) || ai_is_inexact(im))) {
            long bits = numeric_min_inexact_bits(arg);
            mpfr_prec_t out_prec = (bits && bits > 53) ? (mpfr_prec_t)bits : 53;
            Expr* out = ai_eval_complex(re, im, out_prec, false);
            if (out) return out;
        }
    }
#endif /* USE_MPFR */

    return NULL; /* leave symbolic (e.g. AiryAi[2], AiryAi[x]) */
}

/* ------------------------------------------------------------------ */
/* AiryAiPrime[z] = Ai'(z)                                            */
/* Shares the unified core with AiryAi (the core returns Ai').        */
/* ------------------------------------------------------------------ */

static Expr* airyaiprime_one_arg(Expr* arg) {
    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return ai_prime_value_at_zero();   /* Ai'(0) = -1/(3^(1/3) Gamma[1/3]) */
    /* Ai'(z) -> 0 as z -> +Infinity (recessive solution decays). At -Infinity
     * Ai' oscillates with *growing* amplitude (~|z|^(1/4)) and has no limit,
     * so -Infinity is deliberately left unevaluated (unlike AiryAi). */
    if (ai_is_symbol(arg, "Infinity"))       return expr_new_integer(0);
    if (ai_is_symbol(arg, "Indeterminate"))  return expr_new_symbol("Indeterminate");

#ifdef USE_MPFR
    /* 2. Machine real. */
    if (arg->type == EXPR_REAL)
        return ai_eval_real(arg, 53, true);

    /* 3. Arbitrary-precision real. */
    if (arg->type == EXPR_MPFR)
        return ai_eval_real(arg, mpfr_get_prec(arg->data.mpfr), true);

    /* 4. Complex argument (Complex[..] with an inexact part). */
    {
        Expr *re, *im;
        if (is_complex(arg, &re, &im) && (ai_is_inexact(re) || ai_is_inexact(im))) {
            long bits = numeric_min_inexact_bits(arg);
            mpfr_prec_t out_prec = (bits && bits > 53) ? (mpfr_prec_t)bits : 53;
            Expr* out = ai_eval_complex(re, im, out_prec, true);
            if (out) return out;
        }
    }
#endif /* USE_MPFR */

    return NULL; /* leave symbolic (e.g. AiryAiPrime[2], AiryAiPrime[x]) */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

/* Mathematica-compatible argx diagnostic; returns NULL so the evaluator
 * leaves the call unevaluated. */
static Expr* airyai_emit_argx(size_t argc) {
    fprintf(stderr,
            "AiryAi::argx: AiryAi called with %zu argument%s; "
            "1 argument is expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_airyai(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return airyai_emit_argx(argc);
    return airyai_one_arg(res->data.function.args[0]);
}

/* Mathematica-compatible argx diagnostic for AiryAiPrime. */
static Expr* airyaiprime_emit_argx(size_t argc) {
    fprintf(stderr,
            "AiryAiPrime::argx: AiryAiPrime called with %zu argument%s; "
            "1 argument is expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* AiryAiPrime[z] = Ai'(z). Full numeric evaluator (real/complex, machine and
 * arbitrary precision) sharing AiryAi's unified core, plus the exact value at
 * 0 and the +Infinity limit. It also carries the derivative rule
 * (calculus/deriv.c); AiryAiPrime[0] feeds the Taylor series of AiryAi at 0. */
Expr* builtin_airyaiprime(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return airyaiprime_emit_argx(argc);
    return airyaiprime_one_arg(res->data.function.args[0]);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void airyai_init(void) {
    symtab_add_builtin("AiryAi", builtin_airyai);
    symtab_get_def("AiryAi")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED | ATTR_READPROTECTED);

    symtab_add_builtin("AiryAiPrime", builtin_airyaiprime);
    symtab_get_def("AiryAiPrime")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED | ATTR_READPROTECTED);
    /* Docstrings live in info.c (info_init). */
}

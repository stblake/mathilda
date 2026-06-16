/* Mathilda -- BesselJ[n, z], the Bessel function of the first kind J_n(z).
 *
 *   BesselJ[n, z]   solution of  z^2 y'' + z y' + (z^2 - n^2) y = 0  that is
 *                   regular at the origin. For non-integer order it has a
 *                   branch cut along the negative real z axis (from the
 *                   (z/2)^n factor).
 *
 * Evaluation is layered so each kind of argument takes the most accurate and
 * cheapest route:
 *
 *   exact z == 0, integer order    ->  J_0(0) = 1, J_n(0) = 0 (n != 0)
 *   integer order, real z          ->  MPFR-native mpfr_jn (correctly rounded)
 *   any numeric (with an inexact   ->  unified complex-MPFR core (ncpx):
 *     argument), real or complex         power series or asymptotic series
 *   everything else                ->  stays symbolic (return NULL), letting
 *                                      DownValues (half-integer -> elementary)
 *                                      and Series/D fire instead.
 *
 * The unified core `bj_core` evaluates J_nu(z) for arbitrary complex order nu
 * and argument z in the shared complex-MPFR toolkit `ncpx` (numeric_complex.h;
 * the reusable successor to the file-local acx/ecx/gcx of airyai.c/erf.c/
 * gamma.c). It routes between two algorithms on r = |z| and the requested
 * output precision P:
 *
 *   - Power series (small/moderate |z|), valid everywhere on the principal
 *     branch (DLMF 10.2.2):
 *         J_nu(z) = Sum_{k>=0} (-1)^k (z/2)^{nu+2k} / (k! Gamma(nu+k+1)).
 *     Computed by a one-Gamma recurrence: t_0 = (z/2)^nu / Gamma(nu+1),
 *     t_k = t_{k-1} * (-(z/2)^2) / (k (nu+k)). The alternating partial sums
 *     reach magnitude ~ I_nu(|z|) ~ e^{|z|} before cancelling, so the core
 *     adds ~ 2|z|/ln2 guard bits to absorb that exactly.
 *
 *   - Asymptotic series (large |z|, |arg z| < pi), DLMF 10.17.3:
 *         J_nu(z) ~ sqrt(2/(pi z)) [cos(w) A(z) - sin(w) B(z)],
 *         w = z - nu pi/2 - pi/4,
 *         A = Sum (-1)^m a_{2m}/z^{2m}, B = Sum (-1)^m a_{2m+1}/z^{2m+1},
 *     with a_0 = 1, a_k = a_{k-1} (mu - (2k-1)^2)/(8k), mu = 4 nu^2; summed to
 *     the optimal (smallest-term) truncation.
 *
 * Integer order is reduced to non-negative order via J_{-n}(z) = (-1)^n J_n(z)
 * (the power series would otherwise hit Gamma poles); the mpfr_jn fast path
 * handles the sign internally.
 *
 * D[BesselJ[n, z], z] = (BesselJ[n-1, z] - BesselJ[n+1, z]) / 2 lives in
 * calculus/deriv.c; Series at 0 / Infinity and the half-integer -> elementary
 * rewrites live in calculus/series.c and src/internal/bessel.m.
 *
 * Attributes: Listable, NumericFunction, Protected, ReadProtected.
 */
#include "bessel.h"
#include "sym_names.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arithmetic.h"      /* is_rational, is_complex, make_complex */
#include "numeric.h"         /* numeric_min_inexact_bits, get_approx_mpfr */
#include "numeric_complex.h" /* ncpx toolkit */
#include "attr.h"
#include "eval.h"            /* eval_and_free */
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

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool bj_is_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if `e` is inexact, descending one level into Complex[re, im]. */
static bool bj_is_inexact(const Expr* e) {
    if (bj_is_inexact_leaf(e)) return true;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im))
        return bj_is_inexact_leaf(re) || bj_is_inexact_leaf(im);
    return false;
}

/* True if `e` is a concrete real number (no imaginary part). */
static bool bj_is_real_numeric(const Expr* e) {
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) return false;
    return expr_is_numeric_like(e);
}

/* If `e` is an integer-valued order representable as a long (an exact Integer,
 * or a Real/MPFR that is exactly integral), store it in *out and return true.
 * Integer order needs the J_{-n} = (-1)^n J_n reduction (the power series hits
 * Gamma poles otherwise) and unlocks the mpfr_jn fast path. */
static bool bj_order_is_long(const Expr* e, long* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (long)e->data.integer; return true; }
    if (e->type == EXPR_REAL) {
        double v = e->data.real;
        if (isfinite(v) && v == floor(v) && fabs(v) < 9.0e15) { *out = (long)v; return true; }
        return false;
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        if (mpfr_integer_p(e->data.mpfr) && mpfr_fits_slong_p(e->data.mpfr, MPFR_RNDN)) {
            *out = mpfr_get_si(e->data.mpfr, MPFR_RNDN);
            return true;
        }
        return false;
    }
#endif
    return false;
}

#ifdef USE_MPFR
#define BRND MPFR_RNDN

/* Set an already-init2'd real mpfr from an exact-or-real leaf. */
static bool bj_set_real(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, BRND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          BRND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        BRND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          BRND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, BRND);
        mpfr_div_si(out, out, (long)d, BRND);
        return true;
    }
    return false;
}

/* Set an already-init'd ncpx from any numeric (real or Complex) leaf. */
static bool bj_set_ncpx(ncpx* out, const Expr* e) {
    bool inexact;
    return get_approx_mpfr(e, out->re, out->im, &inexact);
}

/* out = Gamma(a) for complex a, at working precision wp.
 *   real a  -> mpfr_gamma (handles negative non-integers; NaN at poles),
 *   else    -> the Gamma builtin (its complex Spouge/Lanczos path). */
static void bj_cgamma(ncpx* out, const ncpx* a, mpfr_prec_t wp) {
    if (mpfr_zero_p(a->im)) {
        mpfr_gamma(out->re, a->re, BRND);
        mpfr_set_zero(out->im, +1);
        return;
    }
    Expr* re = expr_new_mpfr_copy(a->re);
    Expr* im = expr_new_mpfr_copy(a->im);
    Expr* arg = make_complex(re, im);
    Expr* call = expr_new_function(expr_new_symbol(SYM_Gamma), (Expr*[]){ arg }, 1);
    Expr* val = eval_and_free(call);
    bool inexact;
    if (!val || !get_approx_mpfr(val, out->re, out->im, &inexact)) {
        mpfr_set_nan(out->re);
        mpfr_set_nan(out->im);
    }
    (void)wp;
    if (val) expr_free(val);
}

/* ------------------------------------------------------------------ */
/* Power series (small/moderate |z|), valid everywhere.               */
/*   J_nu(z) = Sum_k (-1)^k (z/2)^{nu+2k} / (k! Gamma(nu+k+1))         */
/* Output J is init'd by the caller at precision wp.                   */
/* ------------------------------------------------------------------ */
static void bj_series(ncpx* J, const ncpx* nu, const ncpx* z, mpfr_prec_t wp) {
    ncpx zh, mult, t, g, gamarg, tmp;
    ncpx_init(&zh, wp); ncpx_init(&mult, wp); ncpx_init(&t, wp);
    ncpx_init(&g, wp); ncpx_init(&gamarg, wp); ncpx_init(&tmp, wp);

    /* zh = z/2;  mult = -(z/2)^2  (the per-step multiplier numerator). */
    mpfr_div_2ui(zh.re, z->re, 1, BRND);
    mpfr_div_2ui(zh.im, z->im, 1, BRND);
    ncpx_mul(&mult, &zh, &zh, wp);
    ncpx_neg(&mult, &mult);

    /* t_0 = (z/2)^nu / Gamma(nu+1). */
    ncpx_pow(&t, &zh, nu, wp);                 /* (z/2)^nu */
    ncpx_set(&gamarg, nu);
    mpfr_add_ui(gamarg.re, gamarg.re, 1, BRND);  /* nu+1 */
    bj_cgamma(&g, &gamarg, wp);
    ncpx_div(&t, &t, &g, wp);

    ncpx_set(J, &t);                            /* partial sum = t_0 */

    mpfr_t mag, smag, eps, rh;
    mpfr_inits2(wp, mag, smag, eps, rh, (mpfr_ptr)0);
    ncpx_abs(rh, &zh);                          /* |z/2| */
    double rhd = mpfr_get_d(rh, BRND);
    mpfr_set_ui(eps, 1, BRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 8 ? wp - 4 : 1), BRND);

    unsigned long cap = (unsigned long)(4.0 * rhd) + (unsigned long)wp + 1000;

    for (unsigned long k = 1; k <= cap; k++) {
        /* m_k = mult / (k (nu+k)). */
        ncpx_set(&tmp, nu);
        mpfr_add_ui(tmp.re, tmp.re, k, BRND);   /* nu + k */
        mpfr_mul_ui(tmp.re, tmp.re, k, BRND);   /* k(nu+k) */
        mpfr_mul_ui(tmp.im, tmp.im, k, BRND);
        ncpx_div(&t, &t, &tmp, wp);             /* t *= 1/(k(nu+k)) */
        ncpx_mul(&t, &t, &mult, wp);            /* t *= -(z/2)^2 */
        ncpx_add(J, J, &t);

        /* Past the term peak (k > |z/2|), stop when the term is negligible. */
        if ((double)k > rhd + 2.0) {
            ncpx_abs(mag, &t);
            ncpx_abs(smag, J);
            mpfr_mul(smag, smag, eps, BRND);
            if (mpfr_cmp(mag, smag) < 0) break;
        }
    }

    mpfr_clears(mag, smag, eps, rh, (mpfr_ptr)0);
    ncpx_clear(&zh); ncpx_clear(&mult); ncpx_clear(&t);
    ncpx_clear(&g); ncpx_clear(&gamarg); ncpx_clear(&tmp);
}

/* ------------------------------------------------------------------ */
/* Asymptotic series (large |z|, |arg z| < pi), DLMF 10.17.3.         */
/* Output J is init'd by the caller at precision wp.                   */
/* ------------------------------------------------------------------ */
static void bj_asymp(ncpx* J, const ncpx* nu, const ncpx* z, mpfr_prec_t wp) {
    ncpx mu, omega, cosw, sinw, invz, zpow, a, term, A, B, num, two, pref, t1, t2;
    ncpx_init(&mu, wp); ncpx_init(&omega, wp); ncpx_init(&cosw, wp); ncpx_init(&sinw, wp);
    ncpx_init(&invz, wp); ncpx_init(&zpow, wp); ncpx_init(&a, wp); ncpx_init(&term, wp);
    ncpx_init(&A, wp); ncpx_init(&B, wp); ncpx_init(&num, wp); ncpx_init(&two, wp);
    ncpx_init(&pref, wp); ncpx_init(&t1, wp); ncpx_init(&t2, wp);

    mpfr_t pi, half_pi, quarter_pi;
    mpfr_inits2(wp, pi, half_pi, quarter_pi, (mpfr_ptr)0);
    mpfr_const_pi(pi, BRND);
    mpfr_div_2ui(half_pi, pi, 1, BRND);
    mpfr_div_2ui(quarter_pi, pi, 2, BRND);

    /* mu = 4 nu^2. */
    ncpx_mul(&mu, nu, nu, wp);
    mpfr_mul_ui(mu.re, mu.re, 4, BRND);
    mpfr_mul_ui(mu.im, mu.im, 4, BRND);

    /* omega = z - nu pi/2 - pi/4. */
    ncpx_scale(&t1, nu, half_pi);                /* nu pi/2 */
    ncpx_sub(&omega, z, &t1);
    mpfr_sub(omega.re, omega.re, quarter_pi, BRND);
    ncpx_cos(&cosw, &omega, wp);
    ncpx_sin(&sinw, &omega, wp);

    /* invz = 1/z. */
    ncpx_set_ui(&two, 1);
    ncpx_div(&invz, &two, z, wp);

    /* pref = sqrt(2/(pi z)). */
    ncpx_scale(&t1, z, pi);                       /* pi z */
    ncpx_set_ui(&two, 2);
    ncpx_div(&pref, &two, &t1, wp);               /* 2/(pi z) */
    ncpx_sqrt(&pref, &pref, wp);

    /* k = 0: a_0 = 1, contributes to A with sign +. */
    ncpx_set_ui(&a, 1);
    ncpx_set_ui(&zpow, 1);                         /* z^0 */
    ncpx_set_ui(&A, 1);
    ncpx_set_ui(&B, 0);

    mpfr_t mag, prevmag;
    mpfr_inits2(wp, mag, prevmag, (mpfr_ptr)0);
    mpfr_set_inf(prevmag, 1);

    for (unsigned long k = 1; k < 100000; k++) {
        /* a_k = a_{k-1} (mu - (2k-1)^2) / (8k). */
        unsigned long odd = 2 * k - 1;
        ncpx_set(&num, &mu);
        mpfr_sub_ui(num.re, num.re, odd * odd, BRND);
        ncpx_mul(&a, &a, &num, wp);
        mpfr_div_ui(a.re, a.re, 8 * k, BRND);
        mpfr_div_ui(a.im, a.im, 8 * k, BRND);

        ncpx_mul(&zpow, &zpow, &invz, wp);        /* z^{-k} */
        ncpx_mul(&term, &a, &zpow, wp);           /* a_k / z^k */

        ncpx_abs(mag, &term);
        if (mpfr_cmp(mag, prevmag) >= 0) break;   /* optimal truncation */
        mpfr_set(prevmag, mag, BRND);

        /* sign = (-1)^{floor(k/2)}; even k -> A, odd k -> B. */
        bool negate = ((k / 2) & 1) != 0;
        if ((k & 1) == 0) {
            if (negate) ncpx_sub(&A, &A, &term); else ncpx_add(&A, &A, &term);
        } else {
            if (negate) ncpx_sub(&B, &B, &term); else ncpx_add(&B, &B, &term);
        }
    }

    /* J = pref (cos(w) A - sin(w) B). */
    ncpx_mul(&t1, &cosw, &A, wp);
    ncpx_mul(&t2, &sinw, &B, wp);
    ncpx_sub(&t1, &t1, &t2);
    ncpx_mul(J, &pref, &t1, wp);

    mpfr_clears(mag, prevmag, (mpfr_ptr)0);
    mpfr_clears(pi, half_pi, quarter_pi, (mpfr_ptr)0);
    ncpx_clear(&mu); ncpx_clear(&omega); ncpx_clear(&cosw); ncpx_clear(&sinw);
    ncpx_clear(&invz); ncpx_clear(&zpow); ncpx_clear(&a); ncpx_clear(&term);
    ncpx_clear(&A); ncpx_clear(&B); ncpx_clear(&num); ncpx_clear(&two);
    ncpx_clear(&pref); ncpx_clear(&t1); ncpx_clear(&t2);
}

/* ------------------------------------------------------------------ */
/* Unified core: J_nu(z) at output precision P bits.                  */
/* Picks the power series (small |z|) or asymptotic series (large |z|, */
/* away from the negative-real axis). Output J init'd by the caller.   */
/* ------------------------------------------------------------------ */
static void bj_core(const ncpx* nu, const ncpx* z, ncpx* J, mpfr_prec_t P) {
    mpfr_t rm, am;
    mpfr_init2(rm, (P < 64 ? 64 : P));
    mpfr_init2(am, (P < 64 ? 64 : P));
    ncpx_abs(rm, z);
    ncpx_arg(am, z);
    double rd = mpfr_get_d(rm, BRND);
    double arg = mpfr_get_d(am, BRND);
    mpfr_clear(rm); mpfr_clear(am);

    /* z == 0: J_0(0) = 1, J_nu(0) = 0 for Re(nu) > 0, else undefined. */
    if (rd == 0.0) {
        if (mpfr_zero_p(nu->re) && mpfr_zero_p(nu->im)) ncpx_set_ui(J, 1);
        else if (mpfr_sgn(nu->re) > 0)                  ncpx_set_ui(J, 0);
        else { mpfr_set_nan(J->re); mpfr_set_nan(J->im); }
        return;
    }

    /* Smallest |z| at which the asymptotic series reaches P bits: its
     * optimal-truncation error ~ e^{-2|z|}, so 2|z| > P ln2. */
    double rmin = 0.5 * (double)P * M_LN2;
    if (rmin < 4.0) rmin = 4.0;

    if (rd > 1.3 * rmin && fabs(arg) < 0.9 * M_PI) {
        mpfr_prec_t wp = P + 64;
        ncpx j; ncpx_init(&j, wp);
        bj_asymp(&j, nu, z, wp);
        ncpx_set(J, &j);
        ncpx_clear(&j);
    } else {
        /* Power series. Guard absorbs the ~e^{|z|} (up to e^{2|z/2|})
         * partial-sum cancellation; capped against pathological arguments. */
        long guard = 64 + (long)(2.0 * rd / M_LN2);
        mpfr_prec_t wp = P + (mpfr_prec_t)guard;
        if (wp > P + 300000) wp = P + 300000;
        ncpx j; ncpx_init(&j, wp);
        bj_series(&j, nu, z, wp);
        ncpx_set(J, &j);
        ncpx_clear(&j);
    }
}

/* Build a numeric Expr from a real mpfr value at out_prec bits. */
static Expr* bj_real_result(const mpfr_t v, mpfr_prec_t out_prec) {
    if (out_prec <= 53) {
        double d = mpfr_get_d(v, BRND);
        if (isinf(d) && !mpfr_inf_p(v)) return expr_new_mpfr_copy(v);
        return expr_new_real(d);
    }
    Expr* out = expr_new_mpfr_bits(out_prec);
    mpfr_set(out->data.mpfr, v, BRND);
    return out;
}

/* Build a result from an (re, im) mpfr pair: a real leaf if im rounds to
 * zero, otherwise Complex[..]. Promotes to MPFR on machine-precision
 * overflow so a too-large value is retained rather than lost to Inf. */
static Expr* bj_make_result(const mpfr_t re, const mpfr_t im, mpfr_prec_t out_prec) {
    if (mpfr_zero_p(im)) return bj_real_result(re, out_prec);
    Expr *rr, *ii;
    if (out_prec <= 53) {
        double dr = mpfr_get_d(re, BRND), di = mpfr_get_d(im, BRND);
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
        mpfr_set(rr->data.mpfr, re, BRND);
        mpfr_set(ii->data.mpfr, im, BRND);
    }
    return make_complex(rr, ii);
}

/* Numeric evaluation of J_order(z) at out_prec bits. */
static Expr* bj_eval(const Expr* order, const Expr* z, mpfr_prec_t out_prec) {
    mpfr_prec_t wp = (out_prec < 64 ? 64 : out_prec);

    long n;
    bool int_order = bj_order_is_long(order, &n);

    /* Fast path: integer order, real argument -> mpfr_jn (correctly rounded;
     * handles negative n and negative real z internally). */
    if (int_order && bj_is_real_numeric(z)) {
        mpfr_t zr, out;
        mpfr_init2(zr, wp); mpfr_init2(out, wp);
        Expr* r = NULL;
        if (bj_set_real(zr, z)) {
            long an = n < 0 ? -n : n;
            mpfr_jn(out, an, zr, BRND);
            if (n < 0 && (an & 1)) mpfr_neg(out, out, BRND);
            if (!mpfr_nan_p(out)) r = bj_real_result(out, out_prec);
        }
        mpfr_clear(zr); mpfr_clear(out);
        if (r) return r;
    }

    /* General complex-MPFR core. */
    ncpx nu, zc, J;
    ncpx_init(&nu, wp); ncpx_init(&zc, wp); ncpx_init(&J, wp);
    Expr* out = NULL;
    bool ok;
    if (int_order) {
        /* Reduce to non-negative order; J_{-n} = (-1)^n J_n. */
        long an = n < 0 ? -n : n;
        ncpx_set_ui(&nu, (unsigned long)an);
        ok = bj_set_ncpx(&zc, z);
    } else {
        ok = bj_set_ncpx(&nu, order) && bj_set_ncpx(&zc, z);
    }
    if (ok) {
        bj_core(&nu, &zc, &J, out_prec);
        if (int_order && n < 0 && (n & 1)) ncpx_neg(&J, &J);
        if (!mpfr_nan_p(J.re) && !mpfr_nan_p(J.im) &&
            !mpfr_inf_p(J.re) && !mpfr_inf_p(J.im))
            out = bj_make_result(J.re, J.im, out_prec);
    }
    ncpx_clear(&nu); ncpx_clear(&zc); ncpx_clear(&J);
    return out;
}

/* Governing output precision: min over the inexact arguments, floored at 53. */
static mpfr_prec_t bj_out_prec(const Expr* order, const Expr* z) {
    long bo = numeric_min_inexact_bits(order);
    long bz = numeric_min_inexact_bits(z);
    long b;
    if (bo && bz) b = (bo < bz ? bo : bz);
    else          b = bo ? bo : bz;
    if (b < 53) b = 53;
    return (mpfr_prec_t)b;
}

/* ================================================================== */
/* BesselK -- modified Bessel function of the second kind K_nu(z).    */
/*                                                                    */
/* MPFR has no modified-Bessel routine (only mpfr_jn/mpfr_yn), so the */
/* whole thing is summed in the shared ncpx toolkit. Three kernels,   */
/* picked by |z| and whether the order is an integer:                 */
/*                                                                    */
/*   large |z| (any order)      -> asymptotic series (DLMF 10.40.2)   */
/*   small |z|, non-integer nu  -> connection formula via I_{+-nu}    */
/*                                 (DLMF 10.27.4)                      */
/*   small |z|, integer n       -> logarithmic series (DLMF 10.31.1)  */
/* ================================================================== */

/* Modified Bessel I_mu(z) = Sum_k (z/2)^{mu+2k} / (k! Gamma(mu+k+1)).
 * Same one-Gamma recurrence as bj_series but WITHOUT the (-1)^k sign:
 *   t_0 = (z/2)^mu / Gamma(mu+1),  t_k = t_{k-1} (z/2)^2 / (k (mu+k)).
 * Output I is init'd by the caller at precision wp. */
static void bk_iseries(ncpx* I, const ncpx* mu, const ncpx* z, mpfr_prec_t wp) {
    ncpx zh, mult, t, g, gamarg, tmp;
    ncpx_init(&zh, wp); ncpx_init(&mult, wp); ncpx_init(&t, wp);
    ncpx_init(&g, wp); ncpx_init(&gamarg, wp); ncpx_init(&tmp, wp);

    /* zh = z/2;  mult = (z/2)^2. */
    mpfr_div_2ui(zh.re, z->re, 1, BRND);
    mpfr_div_2ui(zh.im, z->im, 1, BRND);
    ncpx_mul(&mult, &zh, &zh, wp);

    /* t_0 = (z/2)^mu / Gamma(mu+1). For a non-negative integer order use
     * repeated multiplication (exact) instead of the complex pow (log/exp):
     * it keeps a real z's result exactly real, avoiding the spurious tiny
     * imaginary part pow would otherwise introduce for negative real z. */
    if (mpfr_zero_p(mu->im) && mpfr_integer_p(mu->re) &&
        mpfr_sgn(mu->re) >= 0 && mpfr_fits_ulong_p(mu->re, BRND)) {
        unsigned long mn = mpfr_get_ui(mu->re, BRND);
        ncpx_set_ui(&t, 1);
        for (unsigned long i = 0; i < mn; i++) ncpx_mul(&t, &t, &zh, wp);
    } else {
        ncpx_pow(&t, &zh, mu, wp);
    }
    ncpx_set(&gamarg, mu);
    mpfr_add_ui(gamarg.re, gamarg.re, 1, BRND);
    bj_cgamma(&g, &gamarg, wp);
    ncpx_div(&t, &t, &g, wp);

    ncpx_set(I, &t);

    mpfr_t mag, smag, eps, rh;
    mpfr_inits2(wp, mag, smag, eps, rh, (mpfr_ptr)0);
    ncpx_abs(rh, &zh);
    double rhd = mpfr_get_d(rh, BRND);
    mpfr_set_ui(eps, 1, BRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 8 ? wp - 4 : 1), BRND);

    unsigned long cap = (unsigned long)(4.0 * rhd) + (unsigned long)wp + 1000;

    for (unsigned long k = 1; k <= cap; k++) {
        /* t *= (z/2)^2 / (k (mu+k)). */
        ncpx_set(&tmp, mu);
        mpfr_add_ui(tmp.re, tmp.re, k, BRND);
        mpfr_mul_ui(tmp.re, tmp.re, k, BRND);
        mpfr_mul_ui(tmp.im, tmp.im, k, BRND);
        ncpx_div(&t, &t, &tmp, wp);
        ncpx_mul(&t, &t, &mult, wp);
        ncpx_add(I, I, &t);

        if ((double)k > rhd + 2.0) {
            ncpx_abs(mag, &t);
            ncpx_abs(smag, I);
            mpfr_mul(smag, smag, eps, BRND);
            if (mpfr_cmp(mag, smag) < 0) break;
        }
    }

    mpfr_clears(mag, smag, eps, rh, (mpfr_ptr)0);
    ncpx_clear(&zh); ncpx_clear(&mult); ncpx_clear(&t);
    ncpx_clear(&g); ncpx_clear(&gamarg); ncpx_clear(&tmp);
}

/* Asymptotic series (large |z|, |arg z| < 3pi/2), DLMF 10.40.2:
 *   K_nu(z) ~ sqrt(pi/(2z)) e^{-z} Sum_{k>=0} a_k / z^k,
 *   a_0 = 1, a_k = a_{k-1} (4 nu^2 - (2k-1)^2) / (8k).
 * Single sum (no cos/sin split); summed to optimal (smallest-term)
 * truncation. Output K is init'd by the caller at precision wp. */
static void bk_asymp(ncpx* K, const ncpx* nu, const ncpx* z, mpfr_prec_t wp) {
    ncpx mu, invz, zpow, a, term, S, num, pref, negz, ez, halfpi;
    ncpx_init(&mu, wp); ncpx_init(&invz, wp); ncpx_init(&zpow, wp);
    ncpx_init(&a, wp); ncpx_init(&term, wp); ncpx_init(&S, wp);
    ncpx_init(&num, wp); ncpx_init(&pref, wp); ncpx_init(&negz, wp);
    ncpx_init(&ez, wp); ncpx_init(&halfpi, wp);

    mpfr_t pi;
    mpfr_init2(pi, wp);
    mpfr_const_pi(pi, BRND);

    /* mu = 4 nu^2. */
    ncpx_mul(&mu, nu, nu, wp);
    mpfr_mul_ui(mu.re, mu.re, 4, BRND);
    mpfr_mul_ui(mu.im, mu.im, 4, BRND);

    /* invz = 1/z. */
    ncpx_set_ui(&num, 1);
    ncpx_div(&invz, &num, z, wp);

    /* pref = sqrt(pi/(2 z)) * e^{-z}. */
    mpfr_set(halfpi.re, pi, BRND);
    mpfr_div_2ui(halfpi.re, halfpi.re, 1, BRND);   /* pi/2 */
    mpfr_set_zero(halfpi.im, +1);
    ncpx_div(&pref, &halfpi, z, wp);               /* (pi/2)/z = pi/(2z) */
    ncpx_sqrt(&pref, &pref, wp);
    ncpx_neg(&negz, z);
    ncpx_exp(&ez, &negz, wp);
    ncpx_mul(&pref, &pref, &ez, wp);

    /* S = Sum a_k z^{-k}; a_0 = 1, zpow = z^0 = 1. */
    ncpx_set_ui(&a, 1);
    ncpx_set_ui(&zpow, 1);
    ncpx_set_ui(&S, 1);

    mpfr_t mag, prevmag;
    mpfr_inits2(wp, mag, prevmag, (mpfr_ptr)0);
    mpfr_set_inf(prevmag, 1);

    for (unsigned long k = 1; k < 100000; k++) {
        /* a_k = a_{k-1} (mu - (2k-1)^2) / (8k). */
        unsigned long odd = 2 * k - 1;
        ncpx_set(&num, &mu);
        mpfr_sub_ui(num.re, num.re, odd * odd, BRND);
        ncpx_mul(&a, &a, &num, wp);
        mpfr_div_ui(a.re, a.re, 8 * k, BRND);
        mpfr_div_ui(a.im, a.im, 8 * k, BRND);

        ncpx_mul(&zpow, &zpow, &invz, wp);         /* z^{-k} */
        ncpx_mul(&term, &a, &zpow, wp);            /* a_k / z^k */

        ncpx_abs(mag, &term);
        if (mpfr_cmp(mag, prevmag) >= 0) break;    /* optimal truncation */
        mpfr_set(prevmag, mag, BRND);
        ncpx_add(&S, &S, &term);
    }

    ncpx_mul(K, &pref, &S, wp);

    mpfr_clears(mag, prevmag, (mpfr_ptr)0);
    mpfr_clear(pi);
    ncpx_clear(&mu); ncpx_clear(&invz); ncpx_clear(&zpow); ncpx_clear(&a);
    ncpx_clear(&term); ncpx_clear(&S); ncpx_clear(&num); ncpx_clear(&pref);
    ncpx_clear(&negz); ncpx_clear(&ez); ncpx_clear(&halfpi);
}

/* Logarithmic series for integer order n >= 0 (small |z|), DLMF 10.31.1:
 *   K_n(z) = 1/2 (z/2)^{-n} Sum_{k=0}^{n-1} (n-k-1)!/k! (-z^2/4)^k
 *          + (-1)^{n+1} ln(z/2) I_n(z)
 *          + (-1)^n 1/2 (z/2)^n Sum_{k>=0} [psi(k+1)+psi(n+k+1)] /
 *                                          (k!(n+k)!) (z^2/4)^k,
 * with psi(m+1) = -gamma + H_m. Valid (and complex-correct via the
 * (z/2)^{+-n} and ln(z/2) branches) for complex z too. K_{-n} = K_n, so
 * the caller passes |n|. Output K is init'd by the caller at precision wp. */
static void bk_logseries(long n, const ncpx* z, ncpx* K, mpfr_prec_t wp) {
    ncpx zh, zh2, zhn, zhninv, logzh, I_n, t, acc, finite, tmp;
    ncpx_init(&zh, wp); ncpx_init(&zh2, wp); ncpx_init(&zhn, wp);
    ncpx_init(&zhninv, wp); ncpx_init(&logzh, wp); ncpx_init(&I_n, wp);
    ncpx_init(&t, wp); ncpx_init(&acc, wp); ncpx_init(&finite, wp);
    ncpx_init(&tmp, wp);

    /* zh = z/2;  zh2 = (z/2)^2;  zhn = (z/2)^n (exact via repeated mul). */
    mpfr_div_2ui(zh.re, z->re, 1, BRND);
    mpfr_div_2ui(zh.im, z->im, 1, BRND);
    ncpx_mul(&zh2, &zh, &zh, wp);
    ncpx_set_ui(&zhn, 1);
    for (long i = 0; i < n; i++) ncpx_mul(&zhn, &zhn, &zh, wp);
    ncpx_set_ui(&tmp, 1);
    ncpx_div(&zhninv, &tmp, &zhn, wp);             /* (z/2)^{-n} */
    ncpx_log(&logzh, &zh, wp);                     /* ln(z/2) */

    /* gamma. */
    mpfr_t gam, H_k, H_nk, psum, recip;
    mpfr_inits2(wp, gam, H_k, H_nk, psum, recip, (mpfr_ptr)0);
    mpfr_const_euler(gam, BRND);

    /* --- Term 1: 1/2 (z/2)^{-n} Sum_{k=0}^{n-1} (n-k-1)!/k! (-z^2/4)^k. ---
     * term_0 = (n-1)!;  term_k = term_{k-1} (-z^2/4) / (k (n-k)). */
    ncpx_set_ui(&finite, 0);
    if (n > 0) {
        ncpx_set_ui(&t, 1);
        for (long k = 2; k <= n - 1; k++) mpfr_mul_si(t.re, t.re, k, BRND); /* (n-1)! */
        ncpx_set(&finite, &t);
        for (long k = 1; k <= n - 1; k++) {
            ncpx_mul(&t, &t, &zh2, wp);
            ncpx_neg(&t, &t);                       /* * (-z^2/4) */
            mpfr_set_si(recip, (long)(k * (n - k)), BRND);
            mpfr_div(t.re, t.re, recip, BRND);
            mpfr_div(t.im, t.im, recip, BRND);
            ncpx_add(&finite, &finite, &t);
        }
        ncpx_mul(&finite, &finite, &zhninv, wp);    /* * (z/2)^{-n} */
        mpfr_div_2ui(finite.re, finite.re, 1, BRND);/* * 1/2 */
        mpfr_div_2ui(finite.im, finite.im, 1, BRND);
    }

    /* --- Term 2: (-1)^{n+1} ln(z/2) I_n(z). --- */
    ncpx_set_ui(&tmp, (unsigned long)(n < 0 ? -n : n));
    bk_iseries(&I_n, &tmp, z, wp);
    ncpx_mul(&t, &logzh, &I_n, wp);
    if ((n & 1) == 0) ncpx_neg(&t, &t);             /* (-1)^{n+1}: n even -> - */
    ncpx mid; ncpx_init(&mid, wp);
    ncpx_set(&mid, &t);

    /* --- Term 3: (-1)^n 1/2 (z/2)^n Sum_k psum_k (z^2/4)^k/(k!(n+k)!). ---
     * base_0 = 1/n!;  base_k = base_{k-1} (z^2/4)/(k(n+k)).
     * psum_0 = -2gamma + H_n;  psum_k = psum_{k-1} + 1/k + 1/(n+k).
     * H_n = Sum_{j=1}^{n} 1/j. */
    mpfr_set_zero(H_k, +1);                          /* H_0 = 0 */
    mpfr_set_zero(H_nk, +1);
    for (long j = 1; j <= n; j++) {                 /* H_nk = H_n */
        mpfr_set_si(recip, j, BRND);
        mpfr_ui_div(recip, 1, recip, BRND);
        mpfr_add(H_nk, H_nk, recip, BRND);
    }
    mpfr_mul_2ui(psum, gam, 1, BRND);               /* 2 gamma */
    mpfr_neg(psum, psum, BRND);                      /* -2 gamma */
    mpfr_add(psum, psum, H_k, BRND);                 /* + H_0 */
    mpfr_add(psum, psum, H_nk, BRND);                /* + H_n */

    ncpx_set_ui(&t, 1);                              /* base_0 = 1/n! */
    for (long j = 2; j <= n; j++) mpfr_div_si(t.re, t.re, j, BRND);
    ncpx_scale(&acc, &t, psum);                      /* k = 0 contribution */

    mpfr_t mag, smag, eps;
    mpfr_inits2(wp, mag, smag, eps, (mpfr_ptr)0);
    mpfr_set_ui(eps, 1, BRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 8 ? wp - 4 : 1), BRND);
    mpfr_t rh; mpfr_init2(rh, wp);
    ncpx_abs(rh, &zh);
    double rhd = mpfr_get_d(rh, BRND);
    mpfr_clear(rh);

    unsigned long cap = (unsigned long)(4.0 * rhd) + (unsigned long)wp + 1000;
    for (unsigned long k = 1; k <= cap; k++) {
        /* base *= (z^2/4) / (k (n+k)). */
        ncpx_mul(&t, &t, &zh2, wp);
        mpfr_set_si(recip, (long)((long)k * (n + (long)k)), BRND);
        mpfr_div(t.re, t.re, recip, BRND);
        mpfr_div(t.im, t.im, recip, BRND);
        /* psum += 1/k + 1/(n+k). */
        mpfr_set_si(recip, (long)k, BRND); mpfr_ui_div(recip, 1, recip, BRND);
        mpfr_add(psum, psum, recip, BRND);
        mpfr_set_si(recip, n + (long)k, BRND); mpfr_ui_div(recip, 1, recip, BRND);
        mpfr_add(psum, psum, recip, BRND);

        ncpx_scale(&tmp, &t, psum);                  /* contribution_k */
        ncpx_add(&acc, &acc, &tmp);

        if ((double)k > rhd + 2.0) {
            ncpx_abs(mag, &tmp);
            ncpx_abs(smag, &acc);
            mpfr_mul(smag, smag, eps, BRND);
            if (mpfr_cmp(mag, smag) < 0) break;
        }
    }
    ncpx_mul(&acc, &acc, &zhn, wp);                  /* * (z/2)^n */
    mpfr_div_2ui(acc.re, acc.re, 1, BRND);           /* * 1/2 */
    mpfr_div_2ui(acc.im, acc.im, 1, BRND);
    if (n & 1) { ncpx_neg(&acc, &acc); }             /* (-1)^n */

    /* K_n = finite + mid + acc. */
    ncpx_add(K, &finite, &mid);
    ncpx_add(K, K, &acc);

    mpfr_clears(gam, H_k, H_nk, psum, recip, (mpfr_ptr)0);
    mpfr_clears(mag, smag, eps, (mpfr_ptr)0);
    ncpx_clear(&zh); ncpx_clear(&zh2); ncpx_clear(&zhn); ncpx_clear(&zhninv);
    ncpx_clear(&logzh); ncpx_clear(&I_n); ncpx_clear(&t); ncpx_clear(&acc);
    ncpx_clear(&finite); ncpx_clear(&tmp); ncpx_clear(&mid);
}

/* Connection formula for non-integer order (small |z|), DLMF 10.27.4:
 *   K_nu(z) = (pi/2) (I_{-nu}(z) - I_nu(z)) / sin(nu pi).
 * Loses ~ -log2|sin(nu pi)| bits near integer nu; the caller computes at a
 * precision that already absorbs that. Output K is init'd at precision wp. */
static void bk_connection(const ncpx* nu, const ncpx* z, ncpx* K, mpfr_prec_t wp) {
    ncpx ineg, ipos, negnu, nupi, s, diff;
    ncpx_init(&ineg, wp); ncpx_init(&ipos, wp); ncpx_init(&negnu, wp);
    ncpx_init(&nupi, wp); ncpx_init(&s, wp); ncpx_init(&diff, wp);

    mpfr_t pi; mpfr_init2(pi, wp);
    mpfr_const_pi(pi, BRND);

    ncpx_neg(&negnu, nu);
    bk_iseries(&ineg, &negnu, z, wp);     /* I_{-nu} */
    bk_iseries(&ipos, nu,     z, wp);     /* I_{nu}  */
    ncpx_sub(&diff, &ineg, &ipos);

    ncpx_scale(&nupi, nu, pi);            /* nu pi */
    ncpx_sin(&s, &nupi, wp);

    if (mpfr_zero_p(s.re) && mpfr_zero_p(s.im)) {
        mpfr_set_nan(K->re); mpfr_set_nan(K->im);
    } else {
        ncpx_div(K, &diff, &s, wp);       /* (I_{-nu}-I_nu)/sin(nu pi) */
        mpfr_div_2ui(K->re, K->re, 1, BRND);  /* not yet -- need * pi/2 */
        mpfr_div_2ui(K->im, K->im, 1, BRND);
        ncpx_scale(K, K, pi);             /* * pi, with the /2 above => * pi/2 */
    }

    mpfr_clear(pi);
    ncpx_clear(&ineg); ncpx_clear(&ipos); ncpx_clear(&negnu);
    ncpx_clear(&nupi); ncpx_clear(&s); ncpx_clear(&diff);
}

/* Unified core: K_nu(z) at output precision P bits. `int_order` / `an` carry
 * the reduced non-negative integer order (K_{-n}=K_n) when the order is
 * integral. Output K init'd by the caller. */
static void bk_core(const ncpx* nu, const ncpx* z, ncpx* K, mpfr_prec_t P,
                    bool int_order, long an) {
    mpfr_t rm;
    mpfr_init2(rm, (P < 64 ? 64 : P));
    ncpx_abs(rm, z);
    double rd = mpfr_get_d(rm, BRND);
    mpfr_clear(rm);

    /* z == 0 is handled before the core (K diverges); guard anyway. */
    if (rd == 0.0) { mpfr_set_nan(K->re); mpfr_set_nan(K->im); return; }

    /* Asymptotic reaches P bits once its e^{-2|z|}-scale optimal-truncation
     * error is below 2^{-P}; valid on the whole principal sheet (|arg|<3pi/2). */
    double rmin = 0.5 * (double)P * M_LN2;
    if (rmin < 4.0) rmin = 4.0;

    if (rd > 1.3 * rmin) {
        mpfr_prec_t wp = P + 64;
        ncpx k; ncpx_init(&k, wp);
        bk_asymp(&k, nu, z, wp);
        ncpx_set(K, &k);
        ncpx_clear(&k);
        return;
    }

    if (int_order) {
        /* Logarithmic series. Guard absorbs the ~e^{|z|} I_n cancellation. */
        long guard = 64 + (long)(2.0 * rd / M_LN2);
        mpfr_prec_t wp = P + (mpfr_prec_t)guard;
        if (wp > P + 300000) wp = P + 300000;
        ncpx k; ncpx_init(&k, wp);
        bk_logseries(an, z, &k, wp);
        ncpx_set(K, &k);
        ncpx_clear(&k);
        return;
    }

    /* Non-integer: connection formula. Guard absorbs both the I-series
     * e^{|z|} cancellation and the 1/sin(nu pi) loss near integer nu. */
    double sloss = 0.0;
    {
        double a = mpfr_get_d(nu->re, BRND);
        double b = mpfr_get_d(nu->im, BRND);
        double fr = a - floor(a + 0.5);            /* distance to nearest int */
        /* |sin(nu pi)|^2 = sin^2(pi a) cosh^2(pi b) + cos^2(pi a) sinh^2(pi b). */
        double sa = sin(M_PI * fr), ca = cos(M_PI * fr);
        double chb = cosh(M_PI * b), shb = sinh(M_PI * b);
        double smag2 = sa * sa * chb * chb + ca * ca * shb * shb;
        if (smag2 > 0.0 && smag2 < 1.0) sloss = -0.5 * log(smag2) / M_LN2;
    }
    long guard = 64 + (long)(2.0 * rd / M_LN2) + (long)sloss;
    mpfr_prec_t wp = P + (mpfr_prec_t)guard;
    if (wp > P + 300000) wp = P + 300000;
    ncpx k; ncpx_init(&k, wp);
    bk_connection(nu, z, &k, wp);
    ncpx_set(K, &k);
    ncpx_clear(&k);
}

/* Numeric evaluation of K_order(z) at out_prec bits. */
static Expr* bk_eval(const Expr* order, const Expr* z, mpfr_prec_t out_prec) {
    mpfr_prec_t wp = (out_prec < 64 ? 64 : out_prec);

    long n;
    bool int_order = bj_order_is_long(order, &n);
    long an = n < 0 ? -n : n;

    ncpx nu, zc, K;
    ncpx_init(&nu, wp); ncpx_init(&zc, wp); ncpx_init(&K, wp);
    Expr* out = NULL;
    bool ok;
    if (int_order) {
        ncpx_set_ui(&nu, (unsigned long)an);       /* K_{-n} = K_n */
        ok = bj_set_ncpx(&zc, z);
    } else {
        ok = bj_set_ncpx(&nu, order) && bj_set_ncpx(&zc, z);
    }
    if (ok) {
        bk_core(&nu, &zc, &K, out_prec, int_order, an);
        if (!mpfr_nan_p(K.re) && !mpfr_nan_p(K.im) &&
            !mpfr_inf_p(K.re) && !mpfr_inf_p(K.im))
            out = bj_make_result(K.re, K.im, out_prec);
    }
    ncpx_clear(&nu); ncpx_clear(&zc); ncpx_clear(&K);
    return out;
}

/* ================================================================== */
/* BesselI -- modified Bessel function of the first kind I_nu(z).      */
/*                                                                    */
/* I_nu(z) is the solution of z^2 y'' + z y' - (z^2 + n^2) y = 0      */
/* regular at the origin; it grows like e^{|z|} and is even in the     */
/* order for integer index (I_{-n} = I_n). MPFR has no modified-Bessel */
/* routine, so values are summed in the shared ncpx toolkit:           */
/*                                                                    */
/*   small/moderate |z| (any order) -> power series via bk_iseries     */
/*                                     (the modified-Bessel I series    */
/*                                     already used by BesselK)         */
/*   large |z| (any order)          -> asymptotic series (DLMF 10.40.5)*/
/* ================================================================== */

/* Asymptotic series (large |z|, whole principal sheet), DLMF 10.40.5:
 *   I_nu(z) ~ e^z/sqrt(2 pi z) Sum_{k>=0} (-1)^k a_k / z^k
 *           + e^{s(nu+1/2)pi i} e^{-z}/sqrt(2 pi z) Sum_{k>=0} a_k / z^k,
 *   a_0 = 1, a_k = a_{k-1} (4 nu^2 - (2k-1)^2) / (8k);
 * sign s = +1 for arg z > -pi/2, else -1 (covers -pi < arg z <= pi).
 * The second (e^{-z}) term dominates on the left half-plane, so it is
 * required for correctness there, not merely a refinement. The two sums
 * share the a_k recurrence and z^{-k} powers and truncate together at the
 * optimal (smallest-term) point. Output I init'd by the caller at prec wp. */
static void bi_asymp(ncpx* I, const ncpx* nu, const ncpx* z, mpfr_prec_t wp) {
    ncpx mu, invz, zpow, a, term, A, B, num, pref, ez, enz, phase, w, t1, t2;
    ncpx_init(&mu, wp); ncpx_init(&invz, wp); ncpx_init(&zpow, wp);
    ncpx_init(&a, wp); ncpx_init(&term, wp); ncpx_init(&A, wp);
    ncpx_init(&B, wp); ncpx_init(&num, wp); ncpx_init(&pref, wp);
    ncpx_init(&ez, wp); ncpx_init(&enz, wp); ncpx_init(&phase, wp);
    ncpx_init(&w, wp); ncpx_init(&t1, wp); ncpx_init(&t2, wp);

    mpfr_t pi, argz;
    mpfr_inits2(wp, pi, argz, (mpfr_ptr)0);
    mpfr_const_pi(pi, BRND);
    ncpx_arg(argz, z);

    /* mu = 4 nu^2. */
    ncpx_mul(&mu, nu, nu, wp);
    mpfr_mul_ui(mu.re, mu.re, 4, BRND);
    mpfr_mul_ui(mu.im, mu.im, 4, BRND);

    /* invz = 1/z. */
    ncpx_set_ui(&num, 1);
    ncpx_div(&invz, &num, z, wp);

    /* pref = 1/sqrt(2 pi z). */
    ncpx_scale(&t1, z, pi);            /* pi z */
    mpfr_mul_2ui(t1.re, t1.re, 1, BRND);  /* 2 pi z */
    mpfr_mul_2ui(t1.im, t1.im, 1, BRND);
    ncpx_set_ui(&num, 1);
    ncpx_div(&pref, &num, &t1, wp);   /* 1/(2 pi z) */
    ncpx_sqrt(&pref, &pref, wp);

    /* e^{+z}, e^{-z}. */
    ncpx_exp(&ez, z, wp);
    ncpx_neg(&t1, z);
    ncpx_exp(&enz, &t1, wp);

    /* phase = e^{s (nu+1/2) pi i}, s = +1 if arg z > -pi/2 else -1.
     * Let w = nu + 1/2; i*w = (-w.im, w.re); scale by (s pi); exp. */
    {
        double half_pi_d = 0.5 * mpfr_get_d(pi, BRND);
        double s = (mpfr_get_d(argz, BRND) > -half_pi_d) ? 1.0 : -1.0;
        ncpx_set(&w, nu);
        mpfr_add_d(w.re, w.re, 0.5, BRND);     /* nu + 1/2 */
        mpfr_neg(t1.re, w.im, BRND);            /* i*w real part = -w.im */
        mpfr_set(t1.im, w.re, BRND);            /* i*w imag part =  w.re */
        if (s < 0.0) { mpfr_neg(t1.re, t1.re, BRND); mpfr_neg(t1.im, t1.im, BRND); }
        ncpx_scale(&t1, &t1, pi);               /* * pi */
        ncpx_exp(&phase, &t1, wp);
    }

    /* A = Sum (-1)^k a_k z^{-k};  B = Sum a_k z^{-k};  k=0: a_0 = 1. */
    ncpx_set_ui(&a, 1);
    ncpx_set_ui(&zpow, 1);
    ncpx_set_ui(&A, 1);
    ncpx_set_ui(&B, 1);

    mpfr_t mag, prevmag;
    mpfr_inits2(wp, mag, prevmag, (mpfr_ptr)0);
    mpfr_set_inf(prevmag, 1);

    for (unsigned long k = 1; k < 100000; k++) {
        /* a_k = a_{k-1} (mu - (2k-1)^2) / (8k). */
        unsigned long odd = 2 * k - 1;
        ncpx_set(&num, &mu);
        mpfr_sub_ui(num.re, num.re, odd * odd, BRND);
        ncpx_mul(&a, &a, &num, wp);
        mpfr_div_ui(a.re, a.re, 8 * k, BRND);
        mpfr_div_ui(a.im, a.im, 8 * k, BRND);

        ncpx_mul(&zpow, &zpow, &invz, wp);        /* z^{-k} */
        ncpx_mul(&term, &a, &zpow, wp);           /* a_k / z^k */

        ncpx_abs(mag, &term);
        if (mpfr_cmp(mag, prevmag) >= 0) break;   /* optimal truncation */
        mpfr_set(prevmag, mag, BRND);

        ncpx_add(&B, &B, &term);                  /* + a_k/z^k */
        if (k & 1) ncpx_sub(&A, &A, &term);       /* (-1)^k */
        else       ncpx_add(&A, &A, &term);
    }

    /* I = pref (e^z A + phase e^{-z} B). */
    ncpx_mul(&t1, &ez, &A, wp);
    ncpx_mul(&t2, &phase, &enz, wp);
    ncpx_mul(&t2, &t2, &B, wp);
    ncpx_add(&t1, &t1, &t2);
    ncpx_mul(I, &pref, &t1, wp);

    mpfr_clears(mag, prevmag, (mpfr_ptr)0);
    mpfr_clears(pi, argz, (mpfr_ptr)0);
    ncpx_clear(&mu); ncpx_clear(&invz); ncpx_clear(&zpow); ncpx_clear(&a);
    ncpx_clear(&term); ncpx_clear(&A); ncpx_clear(&B); ncpx_clear(&num);
    ncpx_clear(&pref); ncpx_clear(&ez); ncpx_clear(&enz); ncpx_clear(&phase);
    ncpx_clear(&w); ncpx_clear(&t1); ncpx_clear(&t2);
}

/* Unified core: I_nu(z) at output precision P bits. Picks the power series
 * (small/moderate |z|, via the shared bk_iseries) or the asymptotic series
 * (large |z|). Output I init'd by the caller. */
static void bi_core(const ncpx* nu, const ncpx* z, ncpx* I, mpfr_prec_t P) {
    mpfr_t rm;
    mpfr_init2(rm, (P < 64 ? 64 : P));
    ncpx_abs(rm, z);
    double rd = mpfr_get_d(rm, BRND);
    mpfr_clear(rm);

    /* z == 0: I_0(0) = 1, I_nu(0) = 0 for Re(nu) > 0, else undefined. */
    if (rd == 0.0) {
        if (mpfr_zero_p(nu->re) && mpfr_zero_p(nu->im)) ncpx_set_ui(I, 1);
        else if (mpfr_sgn(nu->re) > 0)                  ncpx_set_ui(I, 0);
        else { mpfr_set_nan(I->re); mpfr_set_nan(I->im); }
        return;
    }

    /* Smallest |z| at which the asymptotic reaches P bits: its optimal-
     * truncation error ~ e^{-2|z|} relative to the leading term, so 2|z| > P ln2. */
    double rmin = 0.5 * (double)P * M_LN2;
    if (rmin < 4.0) rmin = 4.0;

    if (rd > 1.3 * rmin) {
        mpfr_prec_t wp = P + 64;
        ncpx i; ncpx_init(&i, wp);
        bi_asymp(&i, nu, z, wp);
        ncpx_set(I, &i);
        ncpx_clear(&i);
    } else {
        /* Power series. Guard absorbs the ~e^{|z|} partial-sum cancellation
         * that can arise for complex z (real positive z does not cancel). */
        long guard = 64 + (long)(2.0 * rd / M_LN2);
        mpfr_prec_t wp = P + (mpfr_prec_t)guard;
        if (wp > P + 300000) wp = P + 300000;
        ncpx i; ncpx_init(&i, wp);
        bk_iseries(&i, nu, z, wp);    /* shared modified-Bessel I power series */
        ncpx_set(I, &i);
        ncpx_clear(&i);
    }
}

/* Numeric evaluation of I_order(z) at out_prec bits. */
static Expr* bi_eval(const Expr* order, const Expr* z, mpfr_prec_t out_prec) {
    mpfr_prec_t wp = (out_prec < 64 ? 64 : out_prec);

    long n;
    bool int_order = bj_order_is_long(order, &n);

    ncpx nu, zc, I;
    ncpx_init(&nu, wp); ncpx_init(&zc, wp); ncpx_init(&I, wp);
    Expr* out = NULL;
    bool ok;
    if (int_order) {
        /* Even in order; I_{-n} = I_n, and the series hits Gamma poles at
         * negative integer order otherwise. */
        long an = n < 0 ? -n : n;
        ncpx_set_ui(&nu, (unsigned long)an);
        ok = bj_set_ncpx(&zc, z);
    } else {
        ok = bj_set_ncpx(&nu, order) && bj_set_ncpx(&zc, z);
    }
    if (ok) {
        /* The asymptotic core carries a complex phase e^{i(nu+1/2)pi}, which
         * seeds sub-precision rounding noise into the imaginary part even for
         * a real argument. The result is mathematically real exactly when the
         * order is real, z is real, AND (the order is an integer, so I_{-n}=I_n
         * is even and stays real on the whole axis, OR z >= 0). A non-integer
         * order with z < 0 lies on the branch cut and IS genuinely complex. */
        bool ord_real = int_order || mpfr_zero_p(nu.im);
        bool z_real   = mpfr_zero_p(zc.im);
        bool force_real = ord_real && z_real && (int_order || mpfr_sgn(zc.re) >= 0);

        bi_core(&nu, &zc, &I, out_prec);
        if (!mpfr_nan_p(I.re) && !mpfr_nan_p(I.im) &&
            !mpfr_inf_p(I.re) && !mpfr_inf_p(I.im)) {
            if (force_real) out = bj_real_result(I.re, out_prec);
            else            out = bj_make_result(I.re, I.im, out_prec);
        }
    }
    ncpx_clear(&nu); ncpx_clear(&zc); ncpx_clear(&I);
    return out;
}

/* ================================================================== */
/* BesselY -- Bessel function of the second kind Y_nu(z).             */
/*                                                                    */
/* Y_nu is the solution of z^2 y'' + z y' + (z^2 - n^2) y = 0 that is */
/* singular (logarithmically, for integer order) at the origin, and   */
/* has a branch cut along the negative real z axis. MPFR exposes      */
/* mpfr_yn for integer order and real z > 0 (used as a fast path);    */
/* everything else is summed in the shared ncpx toolkit. Three        */
/* kernels, picked by |z| and whether the order is an integer:        */
/*                                                                    */
/*   large |z| (any order)      -> asymptotic series (DLMF 10.17.4)   */
/*   small |z|, non-integer nu  -> connection formula via J_{+-nu}    */
/*                                 (DLMF 10.2.3)                       */
/*   small |z|, integer n       -> logarithmic series (DLMF 10.8.1)   */
/*                                                                    */
/* The power-series and asymptotic kernels of BesselJ (bj_series,     */
/* and the bj_asymp shape) are reused directly.                       */
/* ================================================================== */

/* Asymptotic series (large |z|, |arg z| < pi), DLMF 10.17.4:
 *   Y_nu(z) ~ sqrt(2/(pi z)) [sin(w) A(z) + cos(w) B(z)],
 *   w = z - nu pi/2 - pi/4,
 *   A = Sum (-1)^m a_{2m}/z^{2m}, B = Sum (-1)^m a_{2m+1}/z^{2m+1},
 * with a_0 = 1, a_k = a_{k-1} (mu - (2k-1)^2)/(8k), mu = 4 nu^2; summed to
 * the optimal (smallest-term) truncation. This is the companion of bj_asymp:
 * identical A/B series, but the final combination uses (sin A + cos B)
 * instead of J's (cos A - sin B). Output Y init'd by the caller at prec wp. */
static void by_asymp(ncpx* Y, const ncpx* nu, const ncpx* z, mpfr_prec_t wp) {
    ncpx mu, omega, cosw, sinw, invz, zpow, a, term, A, B, num, two, pref, t1, t2;
    ncpx_init(&mu, wp); ncpx_init(&omega, wp); ncpx_init(&cosw, wp); ncpx_init(&sinw, wp);
    ncpx_init(&invz, wp); ncpx_init(&zpow, wp); ncpx_init(&a, wp); ncpx_init(&term, wp);
    ncpx_init(&A, wp); ncpx_init(&B, wp); ncpx_init(&num, wp); ncpx_init(&two, wp);
    ncpx_init(&pref, wp); ncpx_init(&t1, wp); ncpx_init(&t2, wp);

    mpfr_t pi, half_pi, quarter_pi;
    mpfr_inits2(wp, pi, half_pi, quarter_pi, (mpfr_ptr)0);
    mpfr_const_pi(pi, BRND);
    mpfr_div_2ui(half_pi, pi, 1, BRND);
    mpfr_div_2ui(quarter_pi, pi, 2, BRND);

    /* mu = 4 nu^2. */
    ncpx_mul(&mu, nu, nu, wp);
    mpfr_mul_ui(mu.re, mu.re, 4, BRND);
    mpfr_mul_ui(mu.im, mu.im, 4, BRND);

    /* omega = z - nu pi/2 - pi/4. */
    ncpx_scale(&t1, nu, half_pi);                /* nu pi/2 */
    ncpx_sub(&omega, z, &t1);
    mpfr_sub(omega.re, omega.re, quarter_pi, BRND);
    ncpx_cos(&cosw, &omega, wp);
    ncpx_sin(&sinw, &omega, wp);

    /* invz = 1/z. */
    ncpx_set_ui(&two, 1);
    ncpx_div(&invz, &two, z, wp);

    /* pref = sqrt(2/(pi z)). */
    ncpx_scale(&t1, z, pi);                       /* pi z */
    ncpx_set_ui(&two, 2);
    ncpx_div(&pref, &two, &t1, wp);               /* 2/(pi z) */
    ncpx_sqrt(&pref, &pref, wp);

    /* k = 0: a_0 = 1, contributes to A with sign +. */
    ncpx_set_ui(&a, 1);
    ncpx_set_ui(&zpow, 1);                         /* z^0 */
    ncpx_set_ui(&A, 1);
    ncpx_set_ui(&B, 0);

    mpfr_t mag, prevmag;
    mpfr_inits2(wp, mag, prevmag, (mpfr_ptr)0);
    mpfr_set_inf(prevmag, 1);

    for (unsigned long k = 1; k < 100000; k++) {
        /* a_k = a_{k-1} (mu - (2k-1)^2) / (8k). */
        unsigned long odd = 2 * k - 1;
        ncpx_set(&num, &mu);
        mpfr_sub_ui(num.re, num.re, odd * odd, BRND);
        ncpx_mul(&a, &a, &num, wp);
        mpfr_div_ui(a.re, a.re, 8 * k, BRND);
        mpfr_div_ui(a.im, a.im, 8 * k, BRND);

        ncpx_mul(&zpow, &zpow, &invz, wp);        /* z^{-k} */
        ncpx_mul(&term, &a, &zpow, wp);           /* a_k / z^k */

        ncpx_abs(mag, &term);
        if (mpfr_cmp(mag, prevmag) >= 0) break;   /* optimal truncation */
        mpfr_set(prevmag, mag, BRND);

        /* sign = (-1)^{floor(k/2)}; even k -> A, odd k -> B. */
        bool negate = ((k / 2) & 1) != 0;
        if ((k & 1) == 0) {
            if (negate) ncpx_sub(&A, &A, &term); else ncpx_add(&A, &A, &term);
        } else {
            if (negate) ncpx_sub(&B, &B, &term); else ncpx_add(&B, &B, &term);
        }
    }

    /* Y = pref (sin(w) A + cos(w) B). */
    ncpx_mul(&t1, &sinw, &A, wp);
    ncpx_mul(&t2, &cosw, &B, wp);
    ncpx_add(&t1, &t1, &t2);
    ncpx_mul(Y, &pref, &t1, wp);

    mpfr_clears(mag, prevmag, (mpfr_ptr)0);
    mpfr_clears(pi, half_pi, quarter_pi, (mpfr_ptr)0);
    ncpx_clear(&mu); ncpx_clear(&omega); ncpx_clear(&cosw); ncpx_clear(&sinw);
    ncpx_clear(&invz); ncpx_clear(&zpow); ncpx_clear(&a); ncpx_clear(&term);
    ncpx_clear(&A); ncpx_clear(&B); ncpx_clear(&num); ncpx_clear(&two);
    ncpx_clear(&pref); ncpx_clear(&t1); ncpx_clear(&t2);
}

/* Connection formula for non-integer order (small |z|), DLMF 10.2.3:
 *   Y_nu(z) = (J_nu(z) cos(nu pi) - J_{-nu}(z)) / sin(nu pi).
 * Both J's are summed by the shared bj_series power series; loses
 * ~ -log2|sin(nu pi)| bits near integer nu, which the caller's working
 * precision already absorbs. Output Y is init'd by the caller at prec wp. */
static void by_connection(const ncpx* nu, const ncpx* z, ncpx* Y, mpfr_prec_t wp) {
    ncpx jpos, jneg, negnu, nupi, c, s, num;
    ncpx_init(&jpos, wp); ncpx_init(&jneg, wp); ncpx_init(&negnu, wp);
    ncpx_init(&nupi, wp); ncpx_init(&c, wp); ncpx_init(&s, wp); ncpx_init(&num, wp);

    mpfr_t pi; mpfr_init2(pi, wp);
    mpfr_const_pi(pi, BRND);

    bj_series(&jpos, nu, z, wp);          /* J_nu  */
    ncpx_neg(&negnu, nu);
    bj_series(&jneg, &negnu, z, wp);      /* J_{-nu} */

    ncpx_scale(&nupi, nu, pi);            /* nu pi */
    ncpx_cos(&c, &nupi, wp);
    ncpx_sin(&s, &nupi, wp);

    ncpx_mul(&num, &jpos, &c, wp);        /* J_nu cos(nu pi) */
    ncpx_sub(&num, &num, &jneg);          /* - J_{-nu} */

    if (mpfr_zero_p(s.re) && mpfr_zero_p(s.im)) {
        mpfr_set_nan(Y->re); mpfr_set_nan(Y->im);
    } else {
        ncpx_div(Y, &num, &s, wp);
    }

    mpfr_clear(pi);
    ncpx_clear(&jpos); ncpx_clear(&jneg); ncpx_clear(&negnu);
    ncpx_clear(&nupi); ncpx_clear(&c); ncpx_clear(&s); ncpx_clear(&num);
}

/* Logarithmic series for integer order n >= 0 (small |z|), DLMF 10.8.1:
 *   Y_n(z) = -(1/pi) (z/2)^{-n} Sum_{k=0}^{n-1} (n-k-1)!/k! (z^2/4)^k
 *          + (2/pi) ln(z/2) J_n(z)
 *          - (1/pi) (z/2)^n Sum_{k>=0} [psi(k+1)+psi(n+k+1)] /
 *                                      (k!(n+k)!) (-z^2/4)^k,
 * with psi(m+1) = -gamma + H_m. Note vs the BesselK log series (bk_logseries,
 * DLMF 10.31.1): the FINITE sum here carries (+z^2/4)^k (K alternates), the
 * REGULAR sum carries (-z^2/4)^k (K does not), the log term multiplies J_n
 * (not I_n) by 2/pi, and every part scales by -1/pi (or +2/pi). Valid (and
 * complex-correct via the (z/2)^{+-n} and ln(z/2) branches) for complex z too.
 * Y_{-n} = (-1)^n Y_n is applied by the caller; this routine takes |n|. Output
 * Y is init'd by the caller at precision wp. */
static void by_logseries(long n, const ncpx* z, ncpx* Y, mpfr_prec_t wp) {
    ncpx zh, zh2, zhn, zhninv, logzh, J_n, t, acc, finite, tmp;
    ncpx_init(&zh, wp); ncpx_init(&zh2, wp); ncpx_init(&zhn, wp);
    ncpx_init(&zhninv, wp); ncpx_init(&logzh, wp); ncpx_init(&J_n, wp);
    ncpx_init(&t, wp); ncpx_init(&acc, wp); ncpx_init(&finite, wp);
    ncpx_init(&tmp, wp);

    /* zh = z/2;  zh2 = (z/2)^2;  zhn = (z/2)^n (exact via repeated mul). */
    mpfr_div_2ui(zh.re, z->re, 1, BRND);
    mpfr_div_2ui(zh.im, z->im, 1, BRND);
    ncpx_mul(&zh2, &zh, &zh, wp);
    ncpx_set_ui(&zhn, 1);
    for (long i = 0; i < n; i++) ncpx_mul(&zhn, &zhn, &zh, wp);
    ncpx_set_ui(&tmp, 1);
    ncpx_div(&zhninv, &tmp, &zhn, wp);             /* (z/2)^{-n} */
    ncpx_log(&logzh, &zh, wp);                     /* ln(z/2) */

    /* pi, 1/pi, 2/pi, gamma. */
    mpfr_t pi, invpi, twoinvpi, gam, H_k, H_nk, psum, recip;
    mpfr_inits2(wp, pi, invpi, twoinvpi, gam, H_k, H_nk, psum, recip, (mpfr_ptr)0);
    mpfr_const_pi(pi, BRND);
    mpfr_ui_div(invpi, 1, pi, BRND);               /* 1/pi */
    mpfr_mul_2ui(twoinvpi, invpi, 1, BRND);        /* 2/pi */
    mpfr_const_euler(gam, BRND);

    /* --- Term 1: -(1/pi) (z/2)^{-n} Sum_{k=0}^{n-1} (n-k-1)!/k! (z^2/4)^k. ---
     * term_0 = (n-1)!;  term_k = term_{k-1} (z^2/4) / (k (n-k)). */
    ncpx_set_ui(&finite, 0);
    if (n > 0) {
        ncpx_set_ui(&t, 1);
        for (long k = 2; k <= n - 1; k++) mpfr_mul_si(t.re, t.re, k, BRND); /* (n-1)! */
        ncpx_set(&finite, &t);
        for (long k = 1; k <= n - 1; k++) {
            ncpx_mul(&t, &t, &zh2, wp);             /* * (z^2/4)  (no sign flip) */
            mpfr_set_si(recip, (long)(k * (n - k)), BRND);
            mpfr_div(t.re, t.re, recip, BRND);
            mpfr_div(t.im, t.im, recip, BRND);
            ncpx_add(&finite, &finite, &t);
        }
        ncpx_mul(&finite, &finite, &zhninv, wp);    /* * (z/2)^{-n} */
        ncpx_scale(&finite, &finite, invpi);        /* * 1/pi */
        ncpx_neg(&finite, &finite);                 /* * -1 */
    }

    /* --- Term 2: (2/pi) ln(z/2) J_n(z). --- */
    ncpx_set_ui(&tmp, (unsigned long)(n < 0 ? -n : n));
    bj_series(&J_n, &tmp, z, wp);
    ncpx_mul(&t, &logzh, &J_n, wp);
    ncpx_scale(&t, &t, twoinvpi);                   /* * 2/pi */
    ncpx mid; ncpx_init(&mid, wp);
    ncpx_set(&mid, &t);

    /* --- Term 3: -(1/pi) (z/2)^n Sum_k psum_k (-z^2/4)^k/(k!(n+k)!). ---
     * base_0 = 1/n!;  base_k = base_{k-1} (-z^2/4)/(k(n+k)).
     * psum_0 = -2gamma + H_n;  psum_k = psum_{k-1} + 1/k + 1/(n+k).
     * H_n = Sum_{j=1}^{n} 1/j. */
    mpfr_set_zero(H_k, +1);                          /* H_0 = 0 */
    mpfr_set_zero(H_nk, +1);
    for (long j = 1; j <= n; j++) {                 /* H_nk = H_n */
        mpfr_set_si(recip, j, BRND);
        mpfr_ui_div(recip, 1, recip, BRND);
        mpfr_add(H_nk, H_nk, recip, BRND);
    }
    mpfr_mul_2ui(psum, gam, 1, BRND);               /* 2 gamma */
    mpfr_neg(psum, psum, BRND);                      /* -2 gamma */
    mpfr_add(psum, psum, H_k, BRND);                 /* + H_0 */
    mpfr_add(psum, psum, H_nk, BRND);                /* + H_n */

    ncpx_set_ui(&t, 1);                              /* base_0 = 1/n! */
    for (long j = 2; j <= n; j++) mpfr_div_si(t.re, t.re, j, BRND);
    ncpx_scale(&acc, &t, psum);                      /* k = 0 contribution */

    mpfr_t mag, smag, eps;
    mpfr_inits2(wp, mag, smag, eps, (mpfr_ptr)0);
    mpfr_set_ui(eps, 1, BRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 8 ? wp - 4 : 1), BRND);
    mpfr_t rh; mpfr_init2(rh, wp);
    ncpx_abs(rh, &zh);
    double rhd = mpfr_get_d(rh, BRND);
    mpfr_clear(rh);

    unsigned long cap = (unsigned long)(4.0 * rhd) + (unsigned long)wp + 1000;
    for (unsigned long k = 1; k <= cap; k++) {
        /* base *= (-z^2/4) / (k (n+k)). */
        ncpx_mul(&t, &t, &zh2, wp);
        ncpx_neg(&t, &t);                            /* * (-1): (-z^2/4)^k */
        mpfr_set_si(recip, (long)((long)k * (n + (long)k)), BRND);
        mpfr_div(t.re, t.re, recip, BRND);
        mpfr_div(t.im, t.im, recip, BRND);
        /* psum += 1/k + 1/(n+k). */
        mpfr_set_si(recip, (long)k, BRND); mpfr_ui_div(recip, 1, recip, BRND);
        mpfr_add(psum, psum, recip, BRND);
        mpfr_set_si(recip, n + (long)k, BRND); mpfr_ui_div(recip, 1, recip, BRND);
        mpfr_add(psum, psum, recip, BRND);

        ncpx_scale(&tmp, &t, psum);                  /* contribution_k */
        ncpx_add(&acc, &acc, &tmp);

        if ((double)k > rhd + 2.0) {
            ncpx_abs(mag, &tmp);
            ncpx_abs(smag, &acc);
            mpfr_mul(smag, smag, eps, BRND);
            if (mpfr_cmp(mag, smag) < 0) break;
        }
    }
    ncpx_mul(&acc, &acc, &zhn, wp);                  /* * (z/2)^n */
    ncpx_scale(&acc, &acc, invpi);                   /* * 1/pi */
    ncpx_neg(&acc, &acc);                            /* * -1 */

    /* Y_n = finite + mid + acc. */
    ncpx_add(Y, &finite, &mid);
    ncpx_add(Y, Y, &acc);

    mpfr_clears(pi, invpi, twoinvpi, gam, H_k, H_nk, psum, recip, (mpfr_ptr)0);
    mpfr_clears(mag, smag, eps, (mpfr_ptr)0);
    ncpx_clear(&zh); ncpx_clear(&zh2); ncpx_clear(&zhn); ncpx_clear(&zhninv);
    ncpx_clear(&logzh); ncpx_clear(&J_n); ncpx_clear(&t); ncpx_clear(&acc);
    ncpx_clear(&finite); ncpx_clear(&tmp); ncpx_clear(&mid);
}

/* Unified core: Y_nu(z) at output precision P bits. `int_order` / `an` carry
 * the reduced non-negative integer order (the caller applies Y_{-n}=(-1)^n Y_n)
 * when the order is integral. Output Y init'd by the caller. */
static void by_core(const ncpx* nu, const ncpx* z, ncpx* Y, mpfr_prec_t P,
                    bool int_order, long an) {
    mpfr_t rm;
    mpfr_init2(rm, (P < 64 ? 64 : P));
    ncpx_abs(rm, z);
    double rd = mpfr_get_d(rm, BRND);
    mpfr_clear(rm);

    /* z == 0 is handled before the core (Y diverges); guard anyway. */
    if (rd == 0.0) { mpfr_set_nan(Y->re); mpfr_set_nan(Y->im); return; }

    /* Asymptotic reaches P bits once its e^{-2|z|}-scale optimal-truncation
     * error is below 2^{-P}; valid on the principal sheet (|arg z| < pi). */
    double rmin = 0.5 * (double)P * M_LN2;
    if (rmin < 4.0) rmin = 4.0;

    mpfr_t am; mpfr_init2(am, (P < 64 ? 64 : P));
    ncpx_arg(am, z);
    double arg = mpfr_get_d(am, BRND);
    mpfr_clear(am);

    if (rd > 1.3 * rmin && fabs(arg) < 0.9 * M_PI) {
        mpfr_prec_t wp = P + 64;
        ncpx y; ncpx_init(&y, wp);
        by_asymp(&y, nu, z, wp);
        ncpx_set(Y, &y);
        ncpx_clear(&y);
        return;
    }

    if (int_order) {
        /* Logarithmic series. Guard absorbs the ~e^{|z|} J_n cancellation. */
        long guard = 64 + (long)(2.0 * rd / M_LN2);
        mpfr_prec_t wp = P + (mpfr_prec_t)guard;
        if (wp > P + 300000) wp = P + 300000;
        ncpx y; ncpx_init(&y, wp);
        by_logseries(an, z, &y, wp);
        ncpx_set(Y, &y);
        ncpx_clear(&y);
        return;
    }

    /* Non-integer: connection formula. Guard absorbs both the J-series
     * e^{|z|} cancellation and the 1/sin(nu pi) loss near integer nu. */
    double sloss = 0.0;
    {
        double a = mpfr_get_d(nu->re, BRND);
        double b = mpfr_get_d(nu->im, BRND);
        double fr = a - floor(a + 0.5);            /* distance to nearest int */
        /* |sin(nu pi)|^2 = sin^2(pi a) cosh^2(pi b) + cos^2(pi a) sinh^2(pi b). */
        double sa = sin(M_PI * fr), ca = cos(M_PI * fr);
        double chb = cosh(M_PI * b), shb = sinh(M_PI * b);
        double smag2 = sa * sa * chb * chb + ca * ca * shb * shb;
        if (smag2 > 0.0 && smag2 < 1.0) sloss = -0.5 * log(smag2) / M_LN2;
    }
    long guard = 64 + (long)(2.0 * rd / M_LN2) + (long)sloss;
    mpfr_prec_t wp = P + (mpfr_prec_t)guard;
    if (wp > P + 300000) wp = P + 300000;
    ncpx y; ncpx_init(&y, wp);
    by_connection(nu, z, &y, wp);
    ncpx_set(Y, &y);
    ncpx_clear(&y);
}

/* Numeric evaluation of Y_order(z) at out_prec bits. */
static Expr* by_eval(const Expr* order, const Expr* z, mpfr_prec_t out_prec) {
    mpfr_prec_t wp = (out_prec < 64 ? 64 : out_prec);

    long n;
    bool int_order = bj_order_is_long(order, &n);

    /* Fast path: integer order, real argument z > 0 -> mpfr_yn (correctly
     * rounded). Real z <= 0 and complex z fall through to the core: Y has a
     * branch cut on the negative real axis (Y_n(z) is genuinely complex there),
     * which mpfr_yn (real-valued) cannot represent. */
    if (int_order && bj_is_real_numeric(z)) {
        mpfr_t zr, out;
        mpfr_init2(zr, wp); mpfr_init2(out, wp);
        Expr* r = NULL;
        if (bj_set_real(zr, z) && mpfr_sgn(zr) > 0) {
            long an = n < 0 ? -n : n;
            mpfr_yn(out, an, zr, BRND);
            if (n < 0 && (an & 1)) mpfr_neg(out, out, BRND);
            if (!mpfr_nan_p(out)) r = bj_real_result(out, out_prec);
        }
        mpfr_clear(zr); mpfr_clear(out);
        if (r) return r;
    }

    long an = n < 0 ? -n : n;
    ncpx nu, zc, Y;
    ncpx_init(&nu, wp); ncpx_init(&zc, wp); ncpx_init(&Y, wp);
    Expr* out = NULL;
    bool ok;
    if (int_order) {
        /* Reduce to non-negative order; Y_{-n} = (-1)^n Y_n. */
        ncpx_set_ui(&nu, (unsigned long)an);
        ok = bj_set_ncpx(&zc, z);
    } else {
        ok = bj_set_ncpx(&nu, order) && bj_set_ncpx(&zc, z);
    }
    if (ok) {
        by_core(&nu, &zc, &Y, out_prec, int_order, an);
        if (int_order && n < 0 && (n & 1)) ncpx_neg(&Y, &Y);
        if (!mpfr_nan_p(Y.re) && !mpfr_nan_p(Y.im) &&
            !mpfr_inf_p(Y.re) && !mpfr_inf_p(Y.im))
            out = bj_make_result(Y.re, Y.im, out_prec);
    }
    ncpx_clear(&nu); ncpx_clear(&zc); ncpx_clear(&Y);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Dispatch                                                           */
/* ------------------------------------------------------------------ */

/*
 * bessel_z_parity_fold:
 * For INTEGER order n and a superficially-negative argument z, apply the
 * argument-parity reflection  F_n(z) = (-1)^n F_n(-z), where F is BesselJ or
 * BesselI -- both entire in z for integer order. Returns the rewritten
 * expression  (-1)^n F[n, -z]  (with -z canonicalised), or NULL when the fold
 * does not apply (non-integer order, or z not superficially negative).
 *
 * NOT used for BesselK: K_n has a branch cut along the negative real axis, so
 * K_n(-z) carries no clean parity (Mathematica leaves it unevaluated).
 * Non-integer orders are excluded for the same branch-cut reason (J_nu / I_nu
 * are not entire). The negative-integer-order reflection (J_{-n}=(-1)^n J_n,
 * I_{-n}=I_n) is handled separately by the DownValues in src/internal/bessel.m
 * and composes with this fold.
 */
static Expr* bessel_z_parity_fold(const char* head, Expr* order, Expr* z) {
    if (order->type != EXPR_INTEGER) return NULL;
    if (!expr_is_superficially_negative(z)) return NULL;

    /* -z, evaluated so Times[-1, Times[-1, w]] collapses back to w (this also
     * guarantees the recursive call's argument is not superficially negative,
     * so the fold cannot loop). */
    Expr* negz_args[2] = { expr_new_integer(-1), expr_copy(z) };
    Expr* negz = eval_and_free(
        expr_new_function(expr_new_symbol(SYM_Times), negz_args, 2));

    Expr* call_args[2] = { expr_copy(order), negz };
    Expr* call = expr_new_function(expr_new_symbol(head), call_args, 2);

    /* (-1)^n: n even -> +call, n odd -> -call. (C's & on a negative odd long
     * still yields 1 under two's complement, but order here is always the
     * literal integer the user wrote.) */
    if (((long)order->data.integer & 1L) == 0) return call;
    Expr* t_args[2] = { expr_new_integer(-1), call };
    return expr_new_function(expr_new_symbol(SYM_Times), t_args, 2);
}

static Expr* besselj_two_arg(Expr* order, Expr* z) {
    /* Exact special value at the origin (integer order). */
    if (z->type == EXPR_INTEGER && z->data.integer == 0) {
        long n;
        if (order->type == EXPR_INTEGER) {
            n = (long)order->data.integer;
            return expr_new_integer(n == 0 ? 1 : 0);
        }
        return NULL;  /* non-integer order at 0: leave symbolic */
    }

#ifdef USE_MPFR
    /* Numeric path: evaluate when some argument is inexact and both are
     * concrete numbers. All-exact calls stay symbolic so DownValues
     * (half-integer -> elementary) and Series/D can act. */
    if ((bj_is_inexact(order) || bj_is_inexact(z)) &&
        expr_is_numeric_like(order) && expr_is_numeric_like(z)) {
        mpfr_prec_t P = bj_out_prec(order, z);
        Expr* out = bj_eval(order, z, P);
        if (out) return out;
    }
#endif

    /* Argument parity: J_n(-z) = (-1)^n J_n(z) for integer n. Fires only for
     * symbolic / exact-symbolic z (concrete inexact z is handled numerically
     * above). */
    { Expr* f = bessel_z_parity_fold(SYM_BesselJ, order, z); if (f) return f; }

    return NULL;
}

/* Mathematica-compatible argx diagnostic; returns NULL so the evaluator
 * leaves the call unevaluated. */
static Expr* besselj_emit_argx(size_t argc) {
    fprintf(stderr,
            "BesselJ::argrx: BesselJ called with %zu argument%s; "
            "2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_besselj(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return besselj_emit_argx(argc);
    return besselj_two_arg(res->data.function.args[0],
                           res->data.function.args[1]);
}

static Expr* besselk_two_arg(Expr* order, Expr* z) {
    /* Exact pole at the origin: K_0(0) = Infinity (real, directed),
     * K_nu(0) = ComplexInfinity for any other order. */
    if (z->type == EXPR_INTEGER && z->data.integer == 0) {
        if (order->type == EXPR_INTEGER && order->data.integer == 0)
            return expr_new_symbol(SYM_Infinity);
        return expr_new_symbol(SYM_ComplexInfinity);
    }

#ifdef USE_MPFR
    /* Numeric path: evaluate when some argument is inexact and both are
     * concrete numbers. All-exact calls stay symbolic so DownValues
     * (half-integer -> elementary) and Series/D can act. */
    if ((bj_is_inexact(order) || bj_is_inexact(z)) &&
        expr_is_numeric_like(order) && expr_is_numeric_like(z)) {
        mpfr_prec_t P = bj_out_prec(order, z);
        Expr* out = bk_eval(order, z, P);
        if (out) return out;
    }
#endif

    return NULL;
}

static Expr* besselk_emit_argx(size_t argc) {
    fprintf(stderr,
            "BesselK::argrx: BesselK called with %zu argument%s; "
            "2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_besselk(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return besselk_emit_argx(argc);
    return besselk_two_arg(res->data.function.args[0],
                           res->data.function.args[1]);
}

static Expr* besseli_two_arg(Expr* order, Expr* z) {
    /* Exact special value at the origin (integer order): I_0(0) = 1,
     * I_n(0) = 0 for integer n != 0. */
    if (z->type == EXPR_INTEGER && z->data.integer == 0) {
        if (order->type == EXPR_INTEGER) {
            long n = (long)order->data.integer;
            return expr_new_integer(n == 0 ? 1 : 0);
        }
        return NULL;  /* non-integer order at 0: leave symbolic */
    }

#ifdef USE_MPFR
    /* Numeric path: evaluate when some argument is inexact and both are
     * concrete numbers. All-exact calls stay symbolic so DownValues
     * (half-integer -> elementary) and Series/D can act. */
    if ((bj_is_inexact(order) || bj_is_inexact(z)) &&
        expr_is_numeric_like(order) && expr_is_numeric_like(z)) {
        mpfr_prec_t P = bj_out_prec(order, z);
        Expr* out = bi_eval(order, z, P);
        if (out) return out;
    }
#endif

    /* Argument parity: I_n(-z) = (-1)^n I_n(z) for integer n. Fires only for
     * symbolic / exact-symbolic z (concrete inexact z handled numerically). */
    { Expr* f = bessel_z_parity_fold(SYM_BesselI, order, z); if (f) return f; }

    return NULL;
}

static Expr* besseli_emit_argx(size_t argc) {
    fprintf(stderr,
            "BesselI::argrx: BesselI called with %zu argument%s; "
            "2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_besseli(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return besseli_emit_argx(argc);
    return besseli_two_arg(res->data.function.args[0],
                           res->data.function.args[1]);
}

static Expr* bessely_two_arg(Expr* order, Expr* z) {
    /* Exact singularity at the origin (integer order): Y_0(0) = -Infinity
     * (real, directed), Y_n(0) = ComplexInfinity for integer n != 0. A
     * non-integer order is left symbolic (matching besselj_two_arg). */
    if (z->type == EXPR_INTEGER && z->data.integer == 0) {
        if (order->type == EXPR_INTEGER) {
            if (order->data.integer == 0) {
                Expr* args[2] = { expr_new_integer(-1),
                                  expr_new_symbol(SYM_Infinity) };
                return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
            }
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        return NULL;  /* non-integer order at 0: leave symbolic */
    }

#ifdef USE_MPFR
    /* Numeric path: evaluate when some argument is inexact and both are
     * concrete numbers. All-exact calls stay symbolic so DownValues
     * (half-integer -> elementary) and Series/D can act. */
    if ((bj_is_inexact(order) || bj_is_inexact(z)) &&
        expr_is_numeric_like(order) && expr_is_numeric_like(z)) {
        mpfr_prec_t P = bj_out_prec(order, z);
        Expr* out = by_eval(order, z, P);
        if (out) return out;
    }
#endif

    /* No argument-parity fold: Y_n has a branch cut along the negative real
     * axis, so Y_n(-z) carries no clean parity (cf. bessel_z_parity_fold). */
    return NULL;
}

static Expr* bessely_emit_argx(size_t argc) {
    fprintf(stderr,
            "BesselY::argrx: BesselY called with %zu argument%s; "
            "2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_bessely(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return bessely_emit_argx(argc);
    return bessely_two_arg(res->data.function.args[0],
                           res->data.function.args[1]);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void bessel_init(void) {
    symtab_add_builtin("BesselJ", builtin_besselj);
    symtab_get_def("BesselJ")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED | ATTR_READPROTECTED);

    symtab_add_builtin("BesselK", builtin_besselk);
    symtab_get_def("BesselK")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED | ATTR_READPROTECTED);

    symtab_add_builtin("BesselI", builtin_besseli);
    symtab_get_def("BesselI")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED | ATTR_READPROTECTED);

    symtab_add_builtin("BesselY", builtin_bessely);
    symtab_get_def("BesselY")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED | ATTR_READPROTECTED);
    /* Docstrings live in info.c (info_init). */
}

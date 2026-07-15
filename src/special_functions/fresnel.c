/* Mathilda -- the Fresnel integrals (Pi/2-normalized).
 *
 *   FresnelC[z] = Int_0^z Cos[Pi t^2 / 2] dt
 *   FresnelS[z] = Int_0^z Sin[Pi t^2 / 2] dt
 *
 * Both are entire and odd, with no branch cuts. FresnelC and FresnelS share
 * one numeric kernel: the pair (C, S) is computed together and each builtin
 * returns its component. Evaluation is layered so each kind of argument takes
 * the cheapest route:
 *
 *   exact special values     ->  0, +-1/2, +-I/2, Indeterminate
 *   machine / arbitrary real ->  MPFR series (small |x|) or asymptotic (large)
 *   complex (any precision)  ->  the paired A/B ncpx series with guard bits
 *   everything else          ->  stays symbolic (return NULL)
 *
 * Convergent Maclaurin series (valid for all z):
 *
 *   C(z) = Sum_{m>=0} (-1)^m (Pi/2)^(2m)   z^(4m+1) / ((2m)!   (4m+1)),
 *   S(z) = Sum_{m>=0} (-1)^m (Pi/2)^(2m+1) z^(4m+3) / ((2m+1)! (4m+3)).
 *
 * Equivalently A(z) = C(z) + i S(z) = Sum_{k>=0} (i Pi/2)^k z^(2k+1)/(k!(2k+1))
 * and B(z) = C(z) - i S(z) is the same with i -> -i; then C = (A+B)/2 and
 * S = (A-B)/(2i). The real path sums C and S as two real series directly; the
 * complex path sums A and B together in one ncpx loop.
 *
 * The partial sums reach magnitude ~e^((Pi/2)|z|^2) before the O(1)-sized
 * answer emerges, so the MPFR paths add ~(Pi/2)|z|^2/ln2 guard bits to absorb
 * that cancellation exactly. For large real |x| the convergent series is
 * infeasible; there we use the asymptotic expansion (DLMF 7.12)
 *
 *   C(x) = 1/2 + f(x) sin(Pi x^2/2) - g(x) cos(Pi x^2/2),
 *   S(x) = 1/2 - f(x) cos(Pi x^2/2) - g(x) sin(Pi x^2/2),
 *     f(x) ~ (1/(Pi x))   Sum_j (-1)^j (4j-1)!! / (Pi x^2)^(2j),
 *     g(x) ~ (1/(Pi^2 x^3)) Sum_j (-1)^j (4j+1)!! / (Pi x^2)^(2j),
 *
 * summed to the smallest term (optimal truncation). The asymptotic constant
 * 1/2 is the value only in a sector around the real axis (Stokes phenomenon),
 * so it is used for real inputs only; complex inputs always use the convergent
 * series (correct everywhere, merely costlier for large |z|).
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "fresnel.h"
#include "sym_names.h"

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "arithmetic.h"        /* is_rational, make_complex, make_rational, is_complex */
#include "numeric.h"           /* numeric_min_inexact_bits */
#include "attr.h"
#include "eval.h"              /* eval_and_free */
#include "expr.h"
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#include "numeric_complex.h"   /* ncpx toolkit */
#define FRND MPFR_RNDN
#endif

/* M_PI / M_LN2 are POSIX, not C99 -- provide fallbacks (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

static bool fresnel_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

static bool fresnel_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

static bool fresnel_to_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
    int64_t n, d;
    if (is_rational(e, &n, &d)) { *out = (double)n / (double)d; return true; }
    return false;
}

/* True if `e` is -Infinity, represented as Times[-1, Infinity]. */
static bool fresnel_is_neg_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!fresnel_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == -1 && fresnel_is_symbol(b, "Infinity"))
        return true;
    if (b->type == EXPR_INTEGER && b->data.integer == -1 && fresnel_is_symbol(a, "Infinity"))
        return true;
    return false;
}

/* Match Times[Complex[0, s], Infinity] (s = +-1): +-I Infinity. Returns
 * +1 for I Infinity, -1 for -I Infinity, 0 if not of this form. */
static int fresnel_directed_imag_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return 0;
    if (!fresnel_is_symbol(e->data.function.head, "Times")) return 0;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    Expr *cplx = NULL, *inf = NULL;
    if (fresnel_is_symbol(b, "Infinity")) { cplx = a; inf = b; }
    else if (fresnel_is_symbol(a, "Infinity")) { cplx = b; inf = a; }
    if (!inf) return 0;
    Expr *re, *im;
    if (!is_complex(cplx, &re, &im)) return 0;
    if (!(re->type == EXPR_INTEGER && re->data.integer == 0)) return 0;
    if (im->type == EXPR_INTEGER && im->data.integer == 1)  return 1;
    if (im->type == EXPR_INTEGER && im->data.integer == -1) return -1;
    return 0;
}

/* True if `arg` is Times[negative-literal, ...] -- the odd-symmetry trigger. */
static bool fresnel_is_neg_leading_times(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count < 2) return false;
    if (!fresnel_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    if (a->type == EXPR_INTEGER)  return a->data.integer < 0;
    int64_t n, d;
    if (is_rational(a, &n, &d))   return n < 0;
    if (a->type == EXPR_REAL)     return a->data.real < 0.0;
    return false;
}

/* ------------------------------------------------------------------ */
/* Non-MPFR machine fallback (USE_MPFR=0 builds)                      */
/* ------------------------------------------------------------------ */
#ifndef USE_MPFR
/* A(z) = C(z) + i S(z) = Sum_k (i Pi/2)^k z^(2k+1) / (k! (2k+1)). */
static bool fresnel_machine_A(double complex z, double complex* out) {
    double complex ihp = M_PI_2 * I;      /* i Pi/2 */
    double complex z2 = z * z;
    double complex term = z, sum = z;     /* k = 0 */
    double zabs = cabs(z);
    for (int k = 1; k <= 200000; k++) {
        /* term *= (i Pi/2) z^2 (2k-1) / (k (2k+1)) */
        term *= ihp * z2 * (double)(2 * k - 1) / ((double)k * (double)(2 * k + 1));
        sum += term;
        if ((double)k > (M_PI_2 * zabs * zabs) &&
            cabs(term) <= 1e-17 * (cabs(sum) + 1.0)) break;
    }
    if (!isfinite(creal(sum)) || !isfinite(cimag(sum))) return false;
    *out = sum;
    return true;
}

/* B(z) = C(z) - i S(z), the same series with i -> -i. */
static bool fresnel_machine_B(double complex z, double complex* out) {
    double complex ihp = -M_PI_2 * I;     /* -i Pi/2 */
    double complex z2 = z * z;
    double complex term = z, sum = z;
    double zabs = cabs(z);
    for (int k = 1; k <= 200000; k++) {
        term *= ihp * z2 * (double)(2 * k - 1) / ((double)k * (double)(2 * k + 1));
        sum += term;
        if ((double)k > (M_PI_2 * zabs * zabs) &&
            cabs(term) <= 1e-17 * (cabs(sum) + 1.0)) break;
    }
    if (!isfinite(creal(sum)) || !isfinite(cimag(sum))) return false;
    *out = sum;
    return true;
}
#endif /* !USE_MPFR */

#ifdef USE_MPFR
/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool fresnel_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, FRND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          FRND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        FRND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          FRND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, FRND);
        mpfr_div_si(out, out, (long)d, FRND);
        return true;
    }
    return false;
}

/* -------------------- real paths -------------------- */

/* Convergent series for C(x), S(x) into (C, S) at precision wp, for x >= 0.
 *   C: t_0 = x,        t_m = t_{m-1} * (-(Pi/2)^2 x^4 (4m-3)) / ((2m-1)(2m)(4m+1))
 *   S: s_0 = (Pi/2) x^3/3, s_m = s_{m-1} * (-(Pi/2)^2 x^4 (4m-1)) / ((2m)(2m+1)(4m+3))
 * Both summed until past the term peak and decayed below 2^-(wp-8) of it. */
static void fresnel_real_series(mpfr_t C, mpfr_t S, const mpfr_t x,
                                mpfr_prec_t wp, double xabs) {
    mpfr_t tc, ts, x4, K, mag, peak, thr, tmp;
    mpfr_inits2(wp, tc, ts, x4, K, mag, peak, thr, tmp, (mpfr_ptr)0);

    mpfr_const_pi(K, FRND);
    mpfr_div_2ui(K, K, 1, FRND);          /* Pi/2  */
    mpfr_mul(K, K, K, FRND);              /* (Pi/2)^2 */
    mpfr_mul(x4, x, x, FRND);
    mpfr_mul(x4, x4, x4, FRND);           /* x^4 */

    mpfr_set(tc, x, FRND);                /* C term m=0 */
    mpfr_set(C, x, FRND);
    /* s_0 = (Pi/2) x^3 / 3 */
    mpfr_const_pi(ts, FRND);
    mpfr_div_2ui(ts, ts, 1, FRND);        /* Pi/2 */
    mpfr_mul(tmp, x, x, FRND);
    mpfr_mul(tmp, tmp, x, FRND);          /* x^3 */
    mpfr_mul(ts, ts, tmp, FRND);
    mpfr_div_ui(ts, ts, 3, FRND);
    mpfr_set(S, ts, FRND);

    mpfr_abs(peak, tc, FRND);             /* track largest term magnitude */
    mpfr_abs(mag, ts, FRND);
    if (mpfr_cmp(mag, peak) > 0) mpfr_set(peak, mag, FRND);

    long shift = (long)(wp > 8 ? wp - 8 : 1);
    double phase = M_PI_2 * xabs * xabs;
    unsigned long cap = (unsigned long)(2.0 * phase) + (unsigned long)wp + 1000;

    for (unsigned long m = 1; m <= cap; m++) {
        /* tc *= -(K x^4)(4m-3)/((2m-1)(2m)(4m+1)) */
        mpfr_mul(tc, tc, K, FRND);
        mpfr_mul(tc, tc, x4, FRND);
        mpfr_mul_si(tc, tc, (long)(4 * m - 3), FRND);
        mpfr_div_si(tc, tc, (long)(2 * m - 1), FRND);
        mpfr_div_si(tc, tc, (long)(2 * m), FRND);
        mpfr_div_si(tc, tc, (long)(4 * m + 1), FRND);
        mpfr_neg(tc, tc, FRND);
        mpfr_add(C, C, tc, FRND);

        /* ts *= -(K x^4)(4m-1)/((2m)(2m+1)(4m+3)) */
        mpfr_mul(ts, ts, K, FRND);
        mpfr_mul(ts, ts, x4, FRND);
        mpfr_mul_si(ts, ts, (long)(4 * m - 1), FRND);
        mpfr_div_si(ts, ts, (long)(2 * m), FRND);
        mpfr_div_si(ts, ts, (long)(2 * m + 1), FRND);
        mpfr_div_si(ts, ts, (long)(4 * m + 3), FRND);
        mpfr_neg(ts, ts, FRND);
        mpfr_add(S, S, ts, FRND);

        mpfr_abs(mag, tc, FRND);
        if (mpfr_cmp(mag, peak) > 0) mpfr_set(peak, mag, FRND);
        mpfr_abs(tmp, ts, FRND);
        if (mpfr_cmp(tmp, peak) > 0) mpfr_set(peak, tmp, FRND);
        if (mpfr_cmp(tmp, mag) > 0) mpfr_set(mag, tmp, FRND);   /* mag = max term */

        if ((double)m > phase) {
            mpfr_mul_2si(thr, peak, -shift, FRND);
            if (mpfr_cmp(mag, thr) < 0) break;
        }
    }
    mpfr_clears(tc, ts, x4, K, mag, peak, thr, tmp, (mpfr_ptr)0);
}

/* Sum one auxiliary asymptotic series into `acc` (already = t_0 = 1), advancing
 * the running term `t` by -(a)(b)/(w2)^2 each step to its smallest term, where
 * (a,b) = (4j-1, 4j-3) for the f-series or (4j+1, 4j-1) for the g-series. */
static void fresnel_real_aux(mpfr_t acc, mpfr_t t, const mpfr_t w2sq,
                             mpfr_prec_t wp, const mpfr_t thr, int gseries) {
    mpfr_t mag, prev;
    mpfr_inits2(wp, mag, prev, (mpfr_ptr)0);
    mpfr_abs(prev, t, FRND);
    for (unsigned long j = 1; j <= 100000; j++) {
        long a = gseries ? (long)(4 * j + 1) : (long)(4 * j - 1);
        long b = (long)(4 * j - 1);
        if (!gseries) b = (long)(4 * j - 3);
        mpfr_mul_si(t, t, a, FRND);
        mpfr_mul_si(t, t, b, FRND);
        mpfr_div(t, t, w2sq, FRND);
        mpfr_neg(t, t, FRND);
        mpfr_abs(mag, t, FRND);
        if (mpfr_cmp(mag, prev) >= 0) break;       /* terms growing: truncate */
        mpfr_add(acc, acc, t, FRND);
        mpfr_set(prev, mag, FRND);
        if (mpfr_cmp(mag, thr) < 0) break;
    }
    mpfr_clears(mag, prev, (mpfr_ptr)0);
}

/* Asymptotic C(x), S(x) into (C, S) at precision wp, for x > 0. */
static void fresnel_real_asymptotic(mpfr_t C, mpfr_t S, const mpfr_t x,
                                    mpfr_prec_t wp) {
    mpfr_t w2, w2sq, P, R, tP, tR, f, g, phi, sphi, cphi, thr, pi, tmp;
    mpfr_inits2(wp, w2, w2sq, P, R, tP, tR, f, g, phi, sphi, cphi, thr, pi, tmp,
                (mpfr_ptr)0);
    mpfr_const_pi(pi, FRND);

    mpfr_mul(w2, x, x, FRND);
    mpfr_mul(w2, w2, pi, FRND);           /* w2 = Pi x^2 */
    mpfr_mul(w2sq, w2, w2, FRND);         /* (Pi x^2)^2 */

    mpfr_set_ui(thr, 1, FRND);
    mpfr_mul_2si(thr, thr, -(long)wp, FRND);

    mpfr_set_ui(P, 1, FRND);  mpfr_set_ui(tP, 1, FRND);
    fresnel_real_aux(P, tP, w2sq, wp, thr, 0);        /* P = Sum (4j-1)!! ... */
    mpfr_set_ui(R, 1, FRND);  mpfr_set_ui(tR, 1, FRND);
    fresnel_real_aux(R, tR, w2sq, wp, thr, 1);        /* R = Sum (4j+1)!! ... */

    /* f = P/(Pi x), g = R/(Pi^2 x^3) = R/(Pi * w2 * x). */
    mpfr_mul(tmp, pi, x, FRND);           /* Pi x */
    mpfr_div(f, P, tmp, FRND);
    mpfr_mul(tmp, pi, w2, FRND);          /* Pi * (Pi x^2) = Pi^2 x^2 */
    mpfr_mul(tmp, tmp, x, FRND);          /* Pi^2 x^3 */
    mpfr_div(g, R, tmp, FRND);

    /* phi = Pi x^2 / 2 = w2/2 */
    mpfr_div_2ui(phi, w2, 1, FRND);
    mpfr_sin_cos(sphi, cphi, phi, FRND);

    /* C = 1/2 + f sin(phi) - g cos(phi) */
    mpfr_mul(tmp, f, sphi, FRND);
    mpfr_set_d(C, 0.5, FRND);
    mpfr_add(C, C, tmp, FRND);
    mpfr_mul(tmp, g, cphi, FRND);
    mpfr_sub(C, C, tmp, FRND);

    /* S = 1/2 - f cos(phi) - g sin(phi) */
    mpfr_mul(tmp, f, cphi, FRND);
    mpfr_set_d(S, 0.5, FRND);
    mpfr_sub(S, S, tmp, FRND);
    mpfr_mul(tmp, g, sphi, FRND);
    mpfr_sub(S, S, tmp, FRND);

    mpfr_clears(w2, w2sq, P, R, tP, tR, f, g, phi, sphi, cphi, thr, pi, tmp,
                (mpfr_ptr)0);
}

/* Build a real result at out_prec: <= 53 yields a Real, higher yields MPFR. */
static Expr* fresnel_real_result(const mpfr_t v, mpfr_prec_t out_prec) {
    if (out_prec <= 53) return expr_new_real(mpfr_get_d(v, FRND));
    Expr* out = expr_new_mpfr_bits(out_prec);
    mpfr_set(out->data.mpfr, v, FRND);
    return out;
}

/* C(x), S(x) for a numeric real `arg` at out_prec bits, folding via oddness to
 * x >= 0. Writes freshly-allocated Exprs to *outC and *outS. */
static bool fresnel_mpfr_real(const Expr* arg, mpfr_prec_t out_prec,
                              Expr** outC, Expr** outS) {
    mpfr_t x;
    mpfr_init2(x, out_prec + 16);
    if (!fresnel_set_mpfr(x, arg)) { mpfr_clear(x); return false; }
    int sgn = mpfr_sgn(x);
    if (sgn == 0) {
        mpfr_clear(x);
        *outC = out_prec <= 53 ? expr_new_real(0.0) : expr_new_mpfr_bits(out_prec);
        *outS = out_prec <= 53 ? expr_new_real(0.0) : expr_new_mpfr_bits(out_prec);
        return true;
    }
    double xabs = fabs(mpfr_get_d(x, FRND));

    mpfr_t C, S, ax;
    /* Asymptotic is accurate to out_prec once (Pi/2) x^2 >~ out_prec ln2. */
    double asy_x2 = ((double)out_prec + 24.0) * 2.0 * M_LN2 / M_PI + 8.0;
    if (xabs * xabs > asy_x2) {
        mpfr_prec_t wp = out_prec + 64;
        mpfr_inits2(wp, C, S, ax, (mpfr_ptr)0);
        mpfr_abs(ax, x, FRND);
        fresnel_real_asymptotic(C, S, ax, wp);
    } else {
        long guard = 64 + (long)(M_PI_2 * xabs * xabs / M_LN2);
        mpfr_prec_t wp = out_prec + (mpfr_prec_t)guard;
        mpfr_inits2(wp, C, S, ax, (mpfr_ptr)0);
        mpfr_abs(ax, x, FRND);
        fresnel_real_series(C, S, ax, wp, xabs);
    }
    if (sgn < 0) { mpfr_neg(C, C, FRND); mpfr_neg(S, S, FRND); }  /* both odd */

    *outC = fresnel_real_result(C, out_prec);
    *outS = fresnel_real_result(S, out_prec);
    mpfr_clear(C); mpfr_clear(S); mpfr_clear(ax); mpfr_clear(x);
    return true;
}

/* -------------------- complex path (ncpx) -------------------- */

/* out = (i^sign * Pi/2) * in, i.e. multiply by i (sign=+1) or -i (sign=-1)
 * and scale by the real value `pihalf`. Alias-safe. */
static void fresnel_mul_i_scale(ncpx* out, const ncpx* in, const mpfr_t pihalf,
                                int sign) {
    mpfr_t re, im;
    mpfr_init2(re, mpfr_get_prec(out->re));
    mpfr_init2(im, mpfr_get_prec(out->im));
    /* (a + b i) * i = -b + a i ; * (-i) = b - a i */
    if (sign > 0) { mpfr_neg(re, in->im, FRND); mpfr_set(im, in->re, FRND); }
    else          { mpfr_set(re, in->im, FRND); mpfr_neg(im, in->re, FRND); }
    mpfr_mul(out->re, re, pihalf, FRND);
    mpfr_mul(out->im, im, pihalf, FRND);
    mpfr_clear(re); mpfr_clear(im);
}

/* Sum A(w) and B(w) = same series with i -> -i, for arbitrary complex w, into
 * (A, B) at precision wp. term_A_k = term_A_{k-1} (i Pi/2) w^2 (2k-1)/(k(2k+1)).
 * Returns false if it fails to converge within the cap. */
static bool fresnel_ncpx_AB(ncpx* A, ncpx* B, const ncpx* w, mpfr_prec_t wp,
                            double wabs) {
    ncpx w2, tA, tB, scr;
    ncpx_init(&w2, wp); ncpx_init(&tA, wp); ncpx_init(&tB, wp); ncpx_init(&scr, wp);
    ncpx_mul(&w2, w, w, wp);
    ncpx_set(&tA, w); ncpx_set(A, w);      /* k = 0 term = w */
    ncpx_set(&tB, w); ncpx_set(B, w);

    mpfr_t pihalf, mag, peak, thr, m2;
    mpfr_inits2(wp, pihalf, mag, peak, thr, m2, (mpfr_ptr)0);
    mpfr_const_pi(pihalf, FRND);
    mpfr_div_2ui(pihalf, pihalf, 1, FRND);   /* Pi/2 */
    ncpx_abs(peak, w);

    long shift = (long)(wp > 8 ? wp - 8 : 1);
    double phase = M_PI_2 * wabs * wabs;
    unsigned long cap = (unsigned long)(2.0 * phase) + (unsigned long)wp + 1000;
    bool ok = false;
    for (unsigned long k = 1; k <= cap; k++) {
        long num = (long)(2 * k - 1);
        long den1 = (long)k, den2 = (long)(2 * k + 1);
        /* tA = (i Pi/2) w^2 (2k-1)/(k(2k+1)) * tA */
        ncpx_mul(&scr, &tA, &w2, wp);
        mpfr_mul_si(scr.re, scr.re, num, FRND); mpfr_mul_si(scr.im, scr.im, num, FRND);
        mpfr_div_si(scr.re, scr.re, den1, FRND); mpfr_div_si(scr.im, scr.im, den1, FRND);
        mpfr_div_si(scr.re, scr.re, den2, FRND); mpfr_div_si(scr.im, scr.im, den2, FRND);
        fresnel_mul_i_scale(&tA, &scr, pihalf, +1);
        ncpx_add(A, A, &tA);

        ncpx_mul(&scr, &tB, &w2, wp);
        mpfr_mul_si(scr.re, scr.re, num, FRND); mpfr_mul_si(scr.im, scr.im, num, FRND);
        mpfr_div_si(scr.re, scr.re, den1, FRND); mpfr_div_si(scr.im, scr.im, den1, FRND);
        mpfr_div_si(scr.re, scr.re, den2, FRND); mpfr_div_si(scr.im, scr.im, den2, FRND);
        fresnel_mul_i_scale(&tB, &scr, pihalf, -1);
        ncpx_add(B, B, &tB);

        ncpx_abs(mag, &tA);
        ncpx_abs(m2, &tB);
        if (mpfr_cmp(m2, mag) > 0) mpfr_set(mag, m2, FRND);   /* mag = max term */
        if (mpfr_cmp(mag, peak) > 0) mpfr_set(peak, mag, FRND);
        if ((double)k > phase) {
            mpfr_mul_2si(thr, peak, -shift, FRND);
            if (mpfr_cmp(mag, thr) < 0) { ok = true; break; }
        }
    }
    mpfr_clears(pihalf, mag, peak, thr, m2, (mpfr_ptr)0);
    ncpx_clear(&w2); ncpx_clear(&tA); ncpx_clear(&tB); ncpx_clear(&scr);
    return ok;
}

/* Build a complex result from (re, im) at out_prec: <= 53 yields Real parts,
 * higher yields MPFR. make_complex drops a zero imaginary part to a bare real. */
static Expr* fresnel_complex_result(const mpfr_t re, const mpfr_t im,
                                    mpfr_prec_t out_prec) {
    Expr *rr, *ii;
    if (out_prec <= 53) {
        rr = expr_new_real(mpfr_get_d(re, FRND));
        ii = expr_new_real(mpfr_get_d(im, FRND));
    } else {
        rr = expr_new_mpfr_bits(out_prec);
        ii = expr_new_mpfr_bits(out_prec);
        mpfr_set(rr->data.mpfr, re, FRND);
        mpfr_set(ii->data.mpfr, im, FRND);
    }
    return make_complex(rr, ii);
}

/* C(z), S(z) for a numeric complex `arg` at out_prec bits, via the convergent
 * A/B series with guard bits. Folds Re(z) < 0 onto Re >= 0 via oddness. */
static bool fresnel_mpfr_complex(Expr* re, Expr* im, mpfr_prec_t out_prec,
                                 Expr** outC, Expr** outS) {
    double red = 0.0, imd = 0.0;
    (void)fresnel_to_double(re, &red);
    (void)fresnel_to_double(im, &imd);
    int sgn = (red < 0.0) ? -1 : 1;             /* fold to Re >= 0 (C,S odd) */
    double zabs = sqrt(red * red + imd * imd);

    /* Partial sums peak ~e^((Pi/2)|z|^2); add that many guard bits. */
    long guard = 64 + (long)(M_PI_2 * zabs * zabs / M_LN2);
    mpfr_prec_t wp = out_prec + (mpfr_prec_t)guard;

    ncpx w, A, B;
    ncpx_init(&w, wp); ncpx_init(&A, wp); ncpx_init(&B, wp);
    bool done = false;
    if (fresnel_set_mpfr(w.re, re) && fresnel_set_mpfr(w.im, im)) {
        if (sgn < 0) ncpx_neg(&w, &w);
        if (fresnel_ncpx_AB(&A, &B, &w, wp, zabs)) {
            /* C = (A+B)/2 ; S = (A-B)/(2i): S.re = (A-B).im/2, S.im = -(A-B).re/2. */
            mpfr_t cre, cim, sre, sim, dre, dim;
            mpfr_inits2(wp, cre, cim, sre, sim, dre, dim, (mpfr_ptr)0);
            mpfr_add(cre, A.re, B.re, FRND); mpfr_div_2ui(cre, cre, 1, FRND);
            mpfr_add(cim, A.im, B.im, FRND); mpfr_div_2ui(cim, cim, 1, FRND);
            mpfr_sub(dre, A.re, B.re, FRND);
            mpfr_sub(dim, A.im, B.im, FRND);
            mpfr_div_2ui(sre, dim, 1, FRND);            /* (A-B).im / 2 */
            mpfr_div_2ui(sim, dre, 1, FRND);
            mpfr_neg(sim, sim, FRND);                   /* -(A-B).re / 2 */
            if (sgn < 0) {
                mpfr_neg(cre, cre, FRND); mpfr_neg(cim, cim, FRND);
                mpfr_neg(sre, sre, FRND); mpfr_neg(sim, sim, FRND);
            }
            if (!mpfr_nan_p(cre) && !mpfr_nan_p(cim) && !mpfr_inf_p(cre) &&
                !mpfr_inf_p(cim)) {
                *outC = fresnel_complex_result(cre, cim, out_prec);
                *outS = fresnel_complex_result(sre, sim, out_prec);
                done = true;
            }
            mpfr_clears(cre, cim, sre, sim, dre, dim, (mpfr_ptr)0);
        }
    }
    ncpx_clear(&w); ncpx_clear(&A); ncpx_clear(&B);
    return done;
}
#endif /* USE_MPFR */

/* Compute the numeric pair (C, S) for `arg`, or return false to stay symbolic.
 * Callers own the returned Exprs. */
static bool fresnel_numeric(Expr* arg, Expr** outC, Expr** outS) {
    /* Numeric real (machine and arbitrary precision). */
    if (arg->type == EXPR_REAL) {
#ifdef USE_MPFR
        return fresnel_mpfr_real(arg, 53, outC, outS);
#else
        double complex a;
        if (fresnel_machine_A(arg->data.real + 0.0 * I, &a)) {
            *outC = expr_new_real(creal(a));
            *outS = expr_new_real(cimag(a));
            return true;
        }
        return false;
#endif
    }
#ifdef USE_MPFR
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        return fresnel_mpfr_real(arg, prec, outC, outS);
    }
#endif

    /* Numeric complex (Complex[..] with an inexact part). */
    {
        Expr *re, *im;
        if (is_complex(arg, &re, &im) &&
            (fresnel_is_inexact(re) || fresnel_is_inexact(im))) {
#ifdef USE_MPFR
            long bits = numeric_min_inexact_bits(arg);
            mpfr_prec_t out_prec = (bits && bits > 53) ? (mpfr_prec_t)bits : 53;
            return fresnel_mpfr_complex(re, im, out_prec, outC, outS);
#else
            double rr, ii;
            if (fresnel_to_double(re, &rr) && fresnel_to_double(im, &ii)) {
                double complex a, b;
                if (fresnel_machine_A(rr + ii * I, &a) &&
                    fresnel_machine_B(rr + ii * I, &b)) {
                    double complex cc = (a + b) / 2.0;
                    double complex ss = (a - b) / (2.0 * I);
                    *outC = (cimag(cc) == 0.0) ? expr_new_real(creal(cc))
                            : make_complex(expr_new_real(creal(cc)),
                                           expr_new_real(cimag(cc)));
                    *outS = (cimag(ss) == 0.0) ? expr_new_real(creal(ss))
                            : make_complex(expr_new_real(creal(ss)),
                                           expr_new_real(cimag(ss)));
                    return true;
                }
            }
#endif
        }
    }
    return false;
}

/* ------------------------------------------------------------------ */
/* Shared single-argument evaluator (want_S selects FresnelC vs FresnelS) */
/* ------------------------------------------------------------------ */

static Expr* fresnel_one_arg(Expr* arg, int want_S) {
    const char* self = want_S ? SYM_FresnelS : SYM_FresnelC;

    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return expr_new_integer(0);                            /* Fresnel*[0] = 0 */
    if (fresnel_is_symbol(arg, "Infinity"))          return make_rational(1, 2);
    if (fresnel_is_neg_infinity(arg))                return make_rational(-1, 2);
    if (fresnel_is_symbol(arg, "ComplexInfinity"))   return expr_new_symbol(SYM_Indeterminate);
    if (fresnel_is_symbol(arg, "Indeterminate"))     return expr_new_symbol(SYM_Indeterminate);
    {
        int s = fresnel_directed_imag_infinity(arg);          /* +-I Infinity */
        if (s != 0) {
            /* FresnelC[+-I Inf] = +-I/2 ; FresnelS[+-I Inf] = -+I/2. */
            int isign = want_S ? -s : s;
            return make_complex(expr_new_integer(0), make_rational(isign, 2));
        }
    }

    /* 2. Numeric (real or complex). */
    {
        Expr *C = NULL, *S = NULL;
        if (fresnel_numeric(arg, &C, &S)) {
            Expr* keep = want_S ? S : C;
            Expr* drop = want_S ? C : S;
            expr_free(drop);
            return keep;
        }
    }

    /* 3. Symbolic odd symmetry: Fresnel*[-x] -> -Fresnel*[x]. */
    if (fresnel_is_neg_leading_times(arg)) {
        Expr* pos = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1), expr_copy(arg) }, 2));
        Expr* inner = expr_new_function(expr_new_symbol(self), &pos, 1);
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1), inner }, 2));
    }

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry points                                               */
/* ------------------------------------------------------------------ */

static Expr* fresnel_emit_argx(const char* name, size_t argc) {
    fprintf(stderr,
            "%s::argx: %s called with %zu arguments; 1 argument is expected.\n",
            name, name, argc);
    return NULL;
}

Expr* builtin_fresnelc(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 1) return fresnel_one_arg(res->data.function.args[0], 0);
    return fresnel_emit_argx("FresnelC", argc);
}

Expr* builtin_fresnels(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 1) return fresnel_one_arg(res->data.function.args[0], 1);
    return fresnel_emit_argx("FresnelS", argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void fresnel_init(void) {
    symtab_add_builtin("FresnelC", builtin_fresnelc);
    symtab_get_def("FresnelC")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_add_builtin("FresnelS", builtin_fresnels);
    symtab_get_def("FresnelS")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstrings live in info.c (info_init). */
}

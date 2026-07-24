/* Mathilda -- the exponential integral Ei.
 *
 *   ExpIntegralEi[z]   Ei(z) = -PV Int_{-z}^Inf e^-t/t dt
 *
 * Ei has a branch cut along the negative real axis (-Inf, 0); the principal
 * value is taken on the cut. Evaluation is layered so each kind of argument
 * takes the cheapest route:
 *
 *   exact special values    ->  -Infinity, Infinity, 0, +-I Pi, Indeterminate
 *   machine real x > 0       ->  MPFR mpfr_eint (correctly rounded)
 *   machine real x < 0       ->  real convergent series with ln|x|
 *   arbitrary real           ->  the same, at the input precision
 *   complex (any precision)  ->  the convergent series in MPFR with guard bits
 *   everything else          ->  stays symbolic (return NULL)
 *
 * The convergent series (DLMF 6.6.2) is
 *
 *   Ei(z) = gamma + Log(z) + Sum_{k>=1} z^k / (k k!),
 *
 * valid for all z != 0. On the real negative axis the real part ln|x| is used
 * (principal value on the cut); for genuinely complex z the principal Log(z)
 * supplies the +-i Pi jump across the cut. The partial sums of the series can
 * reach magnitude ~e^|z| before the answer (~e^Re z) emerges, so the MPFR path
 * adds (|z| + |Re z|)/ln2 guard bits to absorb that cancellation exactly. For
 * x > 0 MPFR's native mpfr_eint is used directly (fast and correctly rounded,
 * which also covers very-high-precision arguments).
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "expintegralei.h"
#include "sym_names.h"

#include <complex.h>
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

/* M_LN2 is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

/* Euler's gamma to double precision, for the non-MPFR fallback paths. */
#define EI_EULER_GAMMA 0.57721566490153286061

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

/* True if `e` is exactly the symbol `name`. */
static bool ei_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool ei_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* Coerce an exact-or-real leaf to a double. Succeeds for Integer, Real,
 * BigInt, Rational and (for sizing only) MPFR. */
static bool ei_to_double(const Expr* e, double* out) {
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
static bool ei_is_neg_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!ei_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == -1 && ei_is_symbol(b, "Infinity"))
        return true;
    if (b->type == EXPR_INTEGER && b->data.integer == -1 && ei_is_symbol(a, "Infinity"))
        return true;
    return false;
}

/* Match Times[Complex[0, s], Infinity] (s = +-1): +-I Infinity. Returns
 * +1 for I Infinity, -1 for -I Infinity, 0 if not of this form. */
static int ei_directed_imag_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return 0;
    if (!ei_is_symbol(e->data.function.head, "Times")) return 0;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    Expr *cplx = NULL, *inf = NULL;
    if (ei_is_symbol(b, "Infinity")) { cplx = a; inf = b; }
    else if (ei_is_symbol(a, "Infinity")) { cplx = b; inf = a; }
    if (!inf) return 0;
    Expr *re, *im;
    if (!is_complex(cplx, &re, &im)) return 0;
    if (!(re->type == EXPR_INTEGER && re->data.integer == 0)) return 0;
    if (im->type == EXPR_INTEGER && im->data.integer == 1)  return 1;
    if (im->type == EXPR_INTEGER && im->data.integer == -1) return -1;
    return 0;
}

/* Build -Infinity = Times[-1, Infinity]. */
static Expr* ei_make_neg_infinity(void) {
    return expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2);
}

/* Build +-I Pi = Times[Complex[0, s], Pi] for s = +-1. */
static Expr* ei_make_i_pi(int s) {
    Expr* imag = make_complex(expr_new_integer(0), expr_new_integer(s));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ imag, expr_new_symbol(SYM_Pi) }, 2));
}

#ifndef USE_MPFR
/* ------------------------------------------------------------------ */
/* Machine-precision fallback paths (USE_MPFR=0 builds)               */
/* ------------------------------------------------------------------ */

/* Ei(x) for a machine real x != 0. Convergent series for moderate |x|,
 * asymptotic series for large |x|. Returns false if the result is not
 * finite (the caller then leaves the call symbolic). */
static bool ei_machine_real(double x, double* out) {
    double ax = fabs(x);
    double v;
    if (ax <= 45.0) {
        double sum = 0.0, p = 1.0;
        for (int k = 1; k <= 100000; k++) {
            p *= x / k;                  /* p = x^k/k!   */
            double term = p / k;         /* x^k/(k k!)   */
            sum += term;
            if ((double)k > ax && fabs(term) <= 1e-18 * fabs(sum)) break;
        }
        v = EI_EULER_GAMMA + log(ax) + sum;
    } else {
        double s = 1.0, t = 1.0;
        for (int k = 1; k <= (int)ax; k++) {
            t *= (double)k / x;          /* k!/x^k */
            if (fabs(t) <= 1e-18 * fabs(s)) break;
            if (k > 5 && fabs(t) >= 1e3 * fabs(s)) break;  /* asymptotic blow-up */
            s += t;
        }
        v = exp(x) / x * s;
    }
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

/* Ei(z) for a machine double-complex z via the convergent series. */
static bool ei_machine_complex(double complex z, double complex* out) {
    double zabs = cabs(z);
    double complex s = 0.0, p = 1.0;
    for (int k = 1; k <= 200000; k++) {
        p *= z / (double)k;              /* z^k/k!     */
        double complex term = p / (double)k;
        s += term;
        if ((double)k > zabs && cabs(term) <= 1e-17 * cabs(s)) break;
    }
    double complex val = EI_EULER_GAMMA + clog(z) + s;   /* clog: principal */
    if (!isfinite(creal(val)) || !isfinite(cimag(val))) return false;
    *out = val;
    return true;
}
#endif /* !USE_MPFR */

#ifdef USE_MPFR
/* ------------------------------------------------------------------ */
/* Complex-MPFR toolkit (pairs of mpfr_t; no MPC library available).  */
/* Mirrors the `ecx` toolkit in erf.c -- each op runs at an explicit  */
/* working precision and is alias-safe (inputs read into temporaries  */
/* before any output component is written).                           */
/* ------------------------------------------------------------------ */

#define ERND MPFR_RNDN

typedef struct { mpfr_t re, im; } ecx;

static void ecx_init(ecx* z, mpfr_prec_t p) { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void ecx_clear(ecx* z)               { mpfr_clear(z->re);    mpfr_clear(z->im);    }

static void ecx_add(ecx* out, const ecx* a, const ecx* b) {
    mpfr_add(out->re, a->re, b->re, ERND);
    mpfr_add(out->im, a->im, b->im, ERND);
}

/* out = a * b. */
static void ecx_mul(ecx* out, const ecx* a, const ecx* b, mpfr_prec_t p) {
    mpfr_t ac, bd, ad, bc;
    mpfr_inits2(p, ac, bd, ad, bc, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, ERND);
    mpfr_mul(bd, a->im, b->im, ERND);
    mpfr_mul(ad, a->re, b->im, ERND);
    mpfr_mul(bc, a->im, b->re, ERND);
    mpfr_sub(out->re, ac, bd, ERND);
    mpfr_add(out->im, ad, bc, ERND);
    mpfr_clears(ac, bd, ad, bc, (mpfr_ptr)0);
}

/* out = z / r for a real scalar r (in place safe). */
static void ecx_div_d(ecx* out, const ecx* z, double r) {
    mpfr_div_d(out->re, z->re, r, ERND);
    mpfr_div_d(out->im, z->im, r, ERND);
}

/* |z| into mag. */
static void ecx_abs(mpfr_t mag, const ecx* z) { mpfr_hypot(mag, z->re, z->im, ERND); }

/* out = Log(z), principal branch: re = ln|z|, im = atan2(im, re). Alias-safe. */
static void ecx_log(ecx* out, const ecx* z, mpfr_prec_t p) {
    mpfr_t r, ang;
    mpfr_inits2(p, r, ang, (mpfr_ptr)0);
    mpfr_hypot(r, z->re, z->im, ERND);
    mpfr_atan2(ang, z->im, z->re, ERND);
    mpfr_log(out->re, r, ERND);
    mpfr_set(out->im, ang, ERND);
    mpfr_clears(r, ang, (mpfr_ptr)0);
}

/* out = exp(z) = e^Re (cos Im + i sin Im). Alias-safe (reads z into temps). */
static void ecx_exp(ecx* out, const ecx* z, mpfr_prec_t p) {
    mpfr_t er, c, s;
    mpfr_inits2(p, er, c, s, (mpfr_ptr)0);
    mpfr_exp(er, z->re, ERND);
    mpfr_sin_cos(s, c, z->im, ERND);
    mpfr_mul(out->re, er, c, ERND);
    mpfr_mul(out->im, er, s, ERND);
    mpfr_clears(er, c, s, (mpfr_ptr)0);
}

/* out = a / b. Alias-safe (reads inputs into temps before writing out). */
static void ecx_div(ecx* out, const ecx* a, const ecx* b, mpfr_prec_t p) {
    mpfr_t ar, ai, br, bi, den, t;
    mpfr_inits2(p, ar, ai, br, bi, den, t, (mpfr_ptr)0);
    mpfr_set(ar, a->re, ERND); mpfr_set(ai, a->im, ERND);
    mpfr_set(br, b->re, ERND); mpfr_set(bi, b->im, ERND);
    mpfr_mul(den, br, br, ERND);          /* |b|^2 = br^2 + bi^2 */
    mpfr_mul(t,   bi, bi, ERND);
    mpfr_add(den, den, t, ERND);
    mpfr_mul(out->re, ar, br, ERND);      /* re = (ar*br + ai*bi)/|b|^2 */
    mpfr_mul(t,       ai, bi, ERND);
    mpfr_add(out->re, out->re, t, ERND);
    mpfr_div(out->re, out->re, den, ERND);
    mpfr_mul(out->im, ai, br, ERND);      /* im = (ai*br - ar*bi)/|b|^2 */
    mpfr_mul(t,       ar, bi, ERND);
    mpfr_sub(out->im, out->im, t, ERND);
    mpfr_div(out->im, out->im, den, ERND);
    mpfr_clears(ar, ai, br, bi, den, t, (mpfr_ptr)0);
}

/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool ei_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, ERND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          ERND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        ERND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          ERND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, ERND);
        mpfr_div_si(out, out, (long)d, ERND);
        return true;
    }
    return false;
}

/* Real convergent series  Sum_{k>=1} x^k/(k k!) into `sum` (precision wp).
 * Stops once past the term peak (k > |x|) and the term has decayed below
 * 2^-(wp-8) of the largest term seen, so the tail is below the rounding of
 * the result after cancellation. */
static void ei_real_series_sum(mpfr_t sum, const mpfr_t x, mpfr_prec_t wp, double xabs) {
    mpfr_t p, term, mag, peak, thr;
    mpfr_inits2(wp, p, term, mag, peak, thr, (mpfr_ptr)0);
    mpfr_set_zero(sum, 1);
    mpfr_set_ui(p, 1, ERND);            /* x^0/0! */
    mpfr_set_zero(peak, 1);
    long shift = (long)(wp > 8 ? wp - 8 : 1);
    unsigned long cap = (unsigned long)(2.0 * xabs) + (unsigned long)wp + 1000;
    for (unsigned long k = 1; k <= cap; k++) {
        mpfr_mul(p, p, x, ERND);
        mpfr_div_ui(p, p, k, ERND);     /* p = x^k/k!     */
        mpfr_div_ui(term, p, k, ERND);  /* term = x^k/(k k!) */
        mpfr_add(sum, sum, term, ERND);
        mpfr_abs(mag, term, ERND);
        if (mpfr_cmp(mag, peak) > 0) mpfr_set(peak, mag, ERND);
        if ((double)k > xabs) {
            mpfr_mul_2si(thr, peak, -shift, ERND);
            if (mpfr_cmp(mag, thr) < 0) break;
        }
    }
    mpfr_clears(p, term, mag, peak, thr, (mpfr_ptr)0);
}

/* Ei(x) for real x via the divergent asymptotic series (DLMF 6.12.2)
 *
 *   Ei(x) ~ (e^x / x) Sum_{k>=0} k!/x^k,
 *
 * summed up to the smallest term (optimal truncation). Used only when |x| is
 * large enough that the smallest term is below the target precision, so the
 * truncation error is negligible. This is the only feasible route once |x| is
 * so large that the convergent series would demand ~|x|/ln2 guard bits. */
static void ei_real_asymptotic(mpfr_t out, const mpfr_t x, mpfr_prec_t wp, double xabs) {
    mpfr_t term, sum, mag, prev, thr, ex;
    mpfr_inits2(wp, term, sum, mag, prev, thr, ex, (mpfr_ptr)0);
    mpfr_set_ui(term, 1, ERND);             /* k = 0 term = 1 */
    mpfr_set_ui(sum, 1, ERND);
    mpfr_set_ui(prev, 1, ERND);             /* |t_0| = 1 */
    mpfr_set_ui(thr, 1, ERND);
    mpfr_mul_2si(thr, thr, -(long)wp, ERND);  /* 2^-wp */
    unsigned long cap = (unsigned long)xabs + 2;  /* never iterate past the term minimum */
    for (unsigned long k = 1; k <= cap; k++) {
        mpfr_mul_ui(term, term, k, ERND);
        mpfr_div(term, term, x, ERND);      /* term *= k/x */
        mpfr_abs(mag, term, ERND);
        if (mpfr_cmp(mag, prev) >= 0) break; /* terms growing: stop before adding */
        mpfr_add(sum, sum, term, ERND);
        mpfr_set(prev, mag, ERND);
        if (mpfr_cmp(mag, thr) < 0) break;   /* tail below target precision */
    }
    mpfr_exp(ex, x, ERND);
    mpfr_div(out, ex, x, ERND);             /* e^x / x */
    mpfr_mul(out, out, sum, ERND);          /* (e^x/x) * Sum k!/x^k */
    mpfr_clears(term, sum, mag, prev, thr, ex, (mpfr_ptr)0);
}

/* Ei(z) for complex z into already-init'd `out` (precision wp), via the
 * convergent series gamma + Log(z) + Sum_{k>=1} z^k/(k k!). Returns false if
 * the series fails to converge within the cap. */
static bool ei_ecx_series(ecx* out, const ecx* z, mpfr_prec_t wp, double zabs) {
    ecx p, term, s, tmp;
    ecx_init(&p, wp); ecx_init(&term, wp); ecx_init(&s, wp); ecx_init(&tmp, wp);
    mpfr_set_ui(p.re, 1, ERND); mpfr_set_zero(p.im, 1);   /* z^0/0! = 1 */
    mpfr_set_zero(s.re, 1);     mpfr_set_zero(s.im, 1);

    mpfr_t mag, peak, thr;
    mpfr_inits2(wp, mag, peak, thr, (mpfr_ptr)0);
    mpfr_set_zero(peak, 1);
    long shift = (long)(wp > 8 ? wp - 8 : 1);
    unsigned long cap = (unsigned long)(2.0 * zabs) + (unsigned long)wp + 1000;
    bool ok = false;
    for (unsigned long k = 1; k <= cap; k++) {
        ecx_mul(&tmp, &p, z, wp);
        ecx_div_d(&p, &tmp, (double)k);     /* p = z^k/k!        */
        ecx_div_d(&term, &p, (double)k);    /* term = z^k/(k k!) */
        ecx_add(&s, &s, &term);
        ecx_abs(mag, &term);
        if (mpfr_cmp(mag, peak) > 0) mpfr_set(peak, mag, ERND);
        if ((double)k > zabs) {
            mpfr_mul_2si(thr, peak, -shift, ERND);
            if (mpfr_cmp(mag, thr) < 0) { ok = true; break; }
        }
    }

    if (ok) {
        /* out = gamma + Log(z) + s. */
        ecx lg;
        ecx_init(&lg, wp);
        ecx_log(&lg, z, wp);
        mpfr_t euler;
        mpfr_init2(euler, wp);
        mpfr_const_euler(euler, ERND);
        mpfr_add(out->re, euler, lg.re, ERND);
        mpfr_set(out->im, lg.im, ERND);
        ecx_add(out, out, &s);
        mpfr_clear(euler);
        ecx_clear(&lg);
    }

    mpfr_clears(mag, peak, thr, (mpfr_ptr)0);
    ecx_clear(&p); ecx_clear(&term); ecx_clear(&s); ecx_clear(&tmp);
    return ok;
}

/* Ei(z) for complex z into already-init'd `out` (precision wp) via the
 * divergent asymptotic expansion (DLMF 6.12.2 continued off the real axis)
 *
 *   Ei(z) ~ (e^z / z) Sum_{k>=0} k!/z^k  +  i pi sign(Im z),
 *
 * summed to the smallest term. The leading e^z/z is the asymptotic of
 * -E_1(-z); the i pi sign(Im z) constant is the jump of the principal Ei
 * across its branch cut (so Ei(+-I Inf) = +-I Pi). `im_sign` is the sign of
 * Im z. Used only when |z| is large enough that the smallest term is below the
 * target precision -- the regime where the convergent series is infeasible. */
static void ei_ecx_asymptotic(ecx* out, const ecx* z, mpfr_prec_t wp,
                              double zabs, int im_sign) {
    ecx term, sum, scaled, ez, ratio;
    ecx_init(&term, wp); ecx_init(&sum, wp); ecx_init(&scaled, wp);
    ecx_init(&ez, wp);   ecx_init(&ratio, wp);
    mpfr_t mag, prev, thr;
    mpfr_inits2(wp, mag, prev, thr, (mpfr_ptr)0);

    /* sum = Sum_{k>=0} k!/z^k, truncated at the smallest term. */
    mpfr_set_ui(term.re, 1, ERND); mpfr_set_zero(term.im, 1);   /* k = 0 term = 1 */
    mpfr_set_ui(sum.re, 1, ERND);  mpfr_set_zero(sum.im, 1);
    mpfr_set_ui(prev, 1, ERND);                                 /* |t_0| = 1 */
    mpfr_set_ui(thr, 1, ERND);
    mpfr_mul_2si(thr, thr, -(long)wp, ERND);                    /* 2^-wp */
    unsigned long cap = (unsigned long)zabs + 2;  /* never iterate past the minimum */
    for (unsigned long k = 1; k <= cap; k++) {
        mpfr_mul_ui(scaled.re, term.re, k, ERND);  /* scaled = term * k */
        mpfr_mul_ui(scaled.im, term.im, k, ERND);
        ecx_div(&term, &scaled, z, wp);            /* term *= k/z */
        ecx_abs(mag, &term);
        if (mpfr_cmp(mag, prev) >= 0) break;       /* terms growing: stop before adding */
        ecx_add(&sum, &sum, &term);
        mpfr_set(prev, mag, ERND);
        if (mpfr_cmp(mag, thr) < 0) break;         /* tail below target precision */
    }

    ecx_exp(&ez, z, wp);
    ecx_div(&ratio, &ez, z, wp);                   /* e^z / z */
    ecx_mul(out, &ratio, &sum, wp);                /* (e^z/z) * Sum k!/z^k */

    if (im_sign != 0) {                            /* + i pi sign(Im z) */
        mpfr_t pi;
        mpfr_init2(pi, wp);
        mpfr_const_pi(pi, ERND);
        if (im_sign < 0) mpfr_neg(pi, pi, ERND);
        mpfr_add(out->im, out->im, pi, ERND);
        mpfr_clear(pi);
    }

    mpfr_clears(mag, prev, thr, (mpfr_ptr)0);
    ecx_clear(&term); ecx_clear(&sum); ecx_clear(&scaled);
    ecx_clear(&ez);   ecx_clear(&ratio);
}

/* Build a complex result from (re, im) at out_prec: machine precision (<= 53)
 * yields Real parts, higher yields MPFR parts. make_complex drops a zero
 * imaginary part to a bare real. When `force_mpfr` is set the result is always
 * emitted as an MPFR at out_prec, even for out_prec <= 53 -- this preserves an
 * explicit low-precision request (e.g. N[..., 8], ~27 bits) so the value prints
 * the requested number of digits instead of collapsing to a machine double
 * (which prints only ~6 significant figures). */
static Expr* ei_complex_result(const mpfr_t re, const mpfr_t im,
                               mpfr_prec_t out_prec, bool force_mpfr) {
    Expr *rr, *ii;
    if (!force_mpfr && out_prec <= 53) {
        rr = expr_new_real(mpfr_get_d(re, ERND));
        ii = expr_new_real(mpfr_get_d(im, ERND));
    } else {
        rr = expr_new_mpfr_bits(out_prec);
        ii = expr_new_mpfr_bits(out_prec);
        mpfr_set(rr->data.mpfr, re, ERND);
        mpfr_set(ii->data.mpfr, im, ERND);
    }
    return make_complex(rr, ii);
}

/* Ei(x) for a numeric real `arg` at out_prec bits. x > 0 uses mpfr_eint
 * (correctly rounded); x < 0 uses the real convergent series with ln|x| (the
 * principal value on the branch cut). x == 0 yields -Infinity. */
static Expr* ei_mpfr_real(const Expr* arg, mpfr_prec_t out_prec) {
    mpfr_t x;
    mpfr_init2(x, out_prec + 16);
    if (!ei_set_mpfr(x, arg)) { mpfr_clear(x); return NULL; }
    int sgn = mpfr_sgn(x);
    if (sgn == 0) { mpfr_clear(x); return ei_make_neg_infinity(); }

    Expr* out;
    if (sgn > 0) {
        if (out_prec <= 53) {
            mpfr_t r;
            mpfr_init2(r, 53);
            mpfr_eint(r, x, ERND);
            out = expr_new_real(mpfr_get_d(r, ERND));
            mpfr_clear(r);
        } else {
            out = expr_new_mpfr_bits(out_prec);
            mpfr_eint(out->data.mpfr, x, ERND);
        }
    } else {
        double xd = mpfr_get_d(x, ERND);
        double xabs = fabs(xd);
        /* Large |x|: the convergent series would need ~|x|/ln2 guard bits (and
         * overflows the guard arithmetic for astronomically large |x|). Switch
         * to the asymptotic series, which is accurate precisely here -- its
         * smallest term ~ e^-|x| sits below the target precision. */
        if (xabs > (double)(out_prec + 24) * M_LN2 + 8.0) {
            mpfr_prec_t wp = out_prec + 64;
            mpfr_t xv, sum;
            mpfr_inits2(wp, xv, sum, (mpfr_ptr)0);
            ei_set_mpfr(xv, arg);
            ei_real_asymptotic(sum, xv, wp, xabs);
            if (out_prec <= 53) {
                out = expr_new_real(mpfr_get_d(sum, ERND));
            } else {
                out = expr_new_mpfr_bits(out_prec);
                mpfr_set(out->data.mpfr, sum, ERND);
            }
            mpfr_clears(xv, sum, (mpfr_ptr)0);
            mpfr_clear(x);
            return out;
        }
        long guard = 64 + (long)(xabs / M_LN2);
        mpfr_prec_t wp = out_prec + (mpfr_prec_t)guard;
        mpfr_t xv, sum, euler, lnx;
        mpfr_inits2(wp, xv, sum, euler, lnx, (mpfr_ptr)0);
        ei_set_mpfr(xv, arg);
        ei_real_series_sum(sum, xv, wp, xabs);
        mpfr_const_euler(euler, ERND);
        mpfr_abs(lnx, xv, ERND);
        mpfr_log(lnx, lnx, ERND);            /* ln|x| */
        mpfr_add(sum, sum, euler, ERND);
        mpfr_add(sum, sum, lnx, ERND);       /* gamma + ln|x| + series */
        if (out_prec <= 53) {
            out = expr_new_real(mpfr_get_d(sum, ERND));
        } else {
            out = expr_new_mpfr_bits(out_prec);
            mpfr_set(out->data.mpfr, sum, ERND);
        }
        mpfr_clears(xv, sum, euler, lnx, (mpfr_ptr)0);
    }
    mpfr_clear(x);
    return out;
}

/* Ei(z) for a numeric complex `arg` (Complex[..] with at least one inexact
 * part) at out_prec bits, via the MPFR series. Returns the result, or NULL if
 * the parts are not numeric / the series diverges / the result overflows. */
static Expr* ei_mpfr_complex(Expr* re, Expr* im, mpfr_prec_t out_prec,
                             bool force_mpfr) {
    double red = 0.0, imd = 0.0;
    (void)ei_to_double(re, &red);
    (void)ei_to_double(im, &imd);
    double zabs = sqrt(red * red + imd * imd);

    /* Large |z|: the convergent series needs ~|z|/ln2 guard bits and ~2|z|
     * terms -- infeasible, and the guard arithmetic overflows for very large
     * |z|. Use the asymptotic expansion instead; its smallest term ~ e^-|z| is
     * below the target precision exactly in this regime. The threshold also
     * keeps the convergent branch's working precision bounded (no overflow). */
    if (zabs > (double)(out_prec + 24) * M_LN2 + 8.0) {
        mpfr_prec_t wp = out_prec + 64;
        int im_sign = (imd > 0.0) ? 1 : (imd < 0.0 ? -1 : 0);
        ecx z, g;
        ecx_init(&z, wp); ecx_init(&g, wp);
        Expr* out = NULL;
        if (ei_set_mpfr(z.re, re) && ei_set_mpfr(z.im, im)) {
            ei_ecx_asymptotic(&g, &z, wp, zabs, im_sign);
            if (!mpfr_nan_p(g.re) && !mpfr_nan_p(g.im) &&
                !mpfr_inf_p(g.re) && !mpfr_inf_p(g.im))
                out = ei_complex_result(g.re, g.im, out_prec, force_mpfr);
        }
        ecx_clear(&z); ecx_clear(&g);
        return out;
    }

    /* Guard bits absorb the cancellation: terms peak ~e^|z|, the answer is
     * ~e^Re z, so the lost leading digits number (|z| + |Re z|)/ln2 bits. */
    long guard = 64 + (long)((zabs + fabs(red)) / M_LN2);
    mpfr_prec_t wp = out_prec + (mpfr_prec_t)guard;

    ecx z, g;
    ecx_init(&z, wp); ecx_init(&g, wp);
    Expr* out = NULL;
    if (ei_set_mpfr(z.re, re) && ei_set_mpfr(z.im, im) &&
        ei_ecx_series(&g, &z, wp, zabs) &&
        !mpfr_nan_p(g.re) && !mpfr_nan_p(g.im) &&
        !mpfr_inf_p(g.re) && !mpfr_inf_p(g.im)) {
        out = ei_complex_result(g.re, g.im, out_prec, force_mpfr);
    }
    ecx_clear(&z); ecx_clear(&g);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* ExpIntegralEi[z]                                                   */
/* ------------------------------------------------------------------ */

static Expr* ei_one_arg(Expr* arg) {
    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return ei_make_neg_infinity();                       /* Ei[0] = -Inf  */
    if (ei_is_symbol(arg, "Infinity"))         return expr_new_symbol(SYM_Infinity);
    if (ei_is_neg_infinity(arg))               return expr_new_integer(0);  /* Ei[-Inf]=0 */
    if (ei_is_symbol(arg, "ComplexInfinity"))  return expr_new_symbol(SYM_Indeterminate);
    if (ei_is_symbol(arg, "Indeterminate"))    return expr_new_symbol(SYM_Indeterminate);
    {
        int s = ei_directed_imag_infinity(arg);
        if (s != 0) return ei_make_i_pi(s);                  /* Ei[+-I Inf]=+-I Pi */
    }

    /* 2. Numeric real (machine and arbitrary precision). */
    if (arg->type == EXPR_REAL) {
#ifdef USE_MPFR
        return ei_mpfr_real(arg, 53);
#else
        double x;
        if (ei_to_double(arg, &x)) {
            if (x == 0.0) return ei_make_neg_infinity();
            double v;
            if (ei_machine_real(x, &v)) return expr_new_real(v);
        }
        return NULL;
#endif
    }
#ifdef USE_MPFR
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        return ei_mpfr_real(arg, prec);
    }
#endif

    /* 3. Numeric complex (Complex[..] with an inexact part). */
    {
        Expr *re, *im;
        if (is_complex(arg, &re, &im) && (ei_is_inexact(re) || ei_is_inexact(im))) {
#ifdef USE_MPFR
            long bits = numeric_min_inexact_bits(arg);
            /* A genuine arbitrary-precision input (any MPFR leaf) carries an
             * explicit precision request -- e.g. N[Ei[-50 I], 8] numericalises
             * the exact -50 I to a 27-bit MPFR complex. Honour that precision
             * even below 53 bits and emit an MPFR result so the value prints
             * the requested digits; a plain machine-double input (53 bits, no
             * MPFR leaf) still yields a machine-double result. */
            bool force_mpfr = numeric_expr_is_mpfr(arg);
            mpfr_prec_t out_prec;
            if (force_mpfr) {
                out_prec = (mpfr_prec_t)(bits > 0 ? bits : 53);
            } else {
                out_prec = (bits && bits > 53) ? (mpfr_prec_t)bits : 53;
            }
            Expr* out = ei_mpfr_complex(re, im, out_prec, force_mpfr);
            if (out) return out;
#else
            double rr, ii;
            if (ei_to_double(re, &rr) && ei_to_double(im, &ii)) {
                double complex v;
                if (ei_machine_complex(rr + ii * I, &v)) {
                    if (cimag(v) == 0.0) return expr_new_real(creal(v));
                    return make_complex(expr_new_real(creal(v)),
                                        expr_new_real(cimag(v)));
                }
            }
#endif
        }
    }

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

/* Mathematica-compatible argx diagnostic; returns NULL so the evaluator
 * leaves the call unevaluated. */
static Expr* ei_emit_argx(size_t argc) {
    fprintf(stderr,
            "ExpIntegralEi::argx: ExpIntegralEi called with %zu arguments; "
            "1 argument is expected.\n",
            argc);
    return NULL;
}

Expr* builtin_expintegralei(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return ei_one_arg(args[0]);
    return ei_emit_argx(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void expintegralei_init(void) {
    symtab_add_builtin("ExpIntegralEi", builtin_expintegralei);
    symtab_get_def("ExpIntegralEi")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

/* Mathilda -- the imaginary error function.
 *
 *   Erfi[z]   imaginary error function   erfi(z) = erf(i z)/i
 *                                                = (2/sqrt(pi)) Int_0^z e^t^2 dt
 *
 * erfi is an entire function (no branch cuts) and odd in z. Evaluation is
 * layered so each kind of argument takes the cheapest route:
 *
 *   exact special values   ->  0, +-Infinity, +-I (the imaginary-axis limits
 *                              are FINITE: erfi(i y) -> i as y -> +Infinity)
 *   symbolic odd argument   ->  Erfi[-x] = -Erfi[x]
 *   machine / arbitrary real -> the all-positive Maclaurin series
 *                              erfi(x) = (2/sqrt(pi)) Sum x^(2n+1)/(n!(2n+1)),
 *                              evaluated in MPFR. For real x every term shares
 *                              x's sign, so the partial sums climb monotonically
 *                              to the result -- no cancellation, no e^x^2
 *                              prefactor, pure-real arithmetic.
 *   complex (any precision)  -> erfi(z) = -i erf(i z), reusing the
 *                              cancellation-aware erf series (DLMF 7.6.2) in
 *                              MPFR with |z|^2/ln2 guard bits, so even
 *                              machine-precision complex results carry full
 *                              accuracy. Double-precision series is the
 *                              USE_MPFR=0 fallback.
 *   everything else         ->  stays symbolic (return NULL)
 *
 * There is no libm erfi and no mpfr_erfi, so both numeric kernels are
 * hand-rolled here.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "erfi.h"
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

/* M_PI / M_LN2 are POSIX, not C99 -- provide fallbacks (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

/* 2/sqrt(pi), the leading coefficient of erfi. */
#define ERFI_TWO_OVER_SQRT_PI 1.1283791670955126

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

/* True if `e` is exactly the symbol `name`. */
static bool erfi_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool erfi_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* Coerce an exact-or-real leaf to a double (for sizing only). */
static bool erfi_to_double(const Expr* e, double* out) {
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

/* True if `e` is a negative real number literal (Integer/Real/Rational). */
static bool erfi_is_neg_number(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer < 0;
    if (e->type == EXPR_REAL)    return e->data.real < 0.0;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    return mpfr_sgn(e->data.mpfr) < 0;
#endif
    int64_t n, d;
    if (is_rational(e, &n, &d)) return ((n < 0) ^ (d < 0)) != 0;
    return false;
}

/* True if `e` is a Times[c, ...] whose leading factor `c` is a negative
 * real number literal -- the trigger for the odd-symmetry rewrite. */
static bool erfi_is_neg_leading_times(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count < 2) return false;
    if (!erfi_is_symbol(e->data.function.head, "Times")) return false;
    return erfi_is_neg_number(e->data.function.args[0]);
}

/* True if `e` is -Infinity, represented as Times[-1, Infinity]. */
static bool erfi_is_neg_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!erfi_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == -1 && erfi_is_symbol(b, "Infinity"))
        return true;
    if (b->type == EXPR_INTEGER && b->data.integer == -1 && erfi_is_symbol(a, "Infinity"))
        return true;
    return false;
}

/* Match Times[Complex[0, s], Infinity] (s = +-1): +-I Infinity. Returns
 * +1 for I Infinity, -1 for -I Infinity, 0 if not of this form. The erfi
 * limits along the imaginary axis are FINITE: erfi(i y) -> i, erfi(-i y) -> -i
 * as y -> +Infinity, so the caller returns +-I rather than a directed
 * infinity. */
static int erfi_directed_imag_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return 0;
    if (!erfi_is_symbol(e->data.function.head, "Times")) return 0;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    Expr *cplx = NULL, *inf = NULL;
    if (erfi_is_symbol(b, "Infinity")) { cplx = a; inf = b; }
    else if (erfi_is_symbol(a, "Infinity")) { cplx = b; inf = a; }
    if (!inf) return 0;
    Expr *re, *im;
    if (!is_complex(cplx, &re, &im)) return 0;
    if (!(re->type == EXPR_INTEGER && re->data.integer == 0)) return 0;
    if (im->type == EXPR_INTEGER && im->data.integer == 1)  return 1;
    if (im->type == EXPR_INTEGER && im->data.integer == -1) return -1;
    return 0;
}

#ifndef USE_MPFR
/* ------------------------------------------------------------------ */
/* Machine double fallbacks (USE_MPFR=0)                              */
/* ------------------------------------------------------------------ */

/* erfi(x) for a machine double x via the all-positive Maclaurin series. */
static bool erfi_machine_real(double x, double* out) {
    double x2 = x * x;
    double u = x;          /* u_0 */
    double s = x;          /* partial sum */
    for (int n = 1; n < 100000; n++) {
        u *= x2 * (2.0 * n - 1.0) / ((double)n * (2.0 * n + 1.0));
        s += u;
        if ((double)n > x2 && fabs(u) <= 1e-17 * fabs(s)) break;
    }
    double v = ERFI_TWO_OVER_SQRT_PI * s;
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

/* erf(z) for a machine double-complex z via the convergent DLMF 7.6.2 series
 * (mirrors erf.c's fallback). */
static bool erfi_erf_machine_complex(double complex z, double complex* out) {
    double complex z2 = z * z;
    double complex two_z2 = 2.0 * z2;
    double complex t = z;          /* t_0 */
    double complex s = z;          /* partial sum */
    double z2abs = cabs(z2);
    for (int n = 1; n < 100000; n++) {
        t *= two_z2 / (2.0 * n + 1.0);
        s += t;
        if ((double)n > z2abs && cabs(t) <= 1e-17 * cabs(s)) break;
    }
    double complex val = ERFI_TWO_OVER_SQRT_PI * cexp(-z2) * s;
    if (!isfinite(creal(val)) || !isfinite(cimag(val))) return false;
    *out = val;
    return true;
}

/* erfi(z) = -i erf(i z) for a machine double-complex z. */
static bool erfi_machine_complex(double complex z, double complex* out) {
    double complex e;
    if (!erfi_erf_machine_complex(z * I, &e)) return false;
    *out = -I * e;
    return true;
}
#endif /* !USE_MPFR */

#ifdef USE_MPFR
/* ------------------------------------------------------------------ */
/* Complex-MPFR toolkit (pairs of mpfr_t; no MPC library available).  */
/* Mirrors the `ecx` toolkit in erf.c -- each op runs at an explicit   */
/* working precision and is alias-safe.                               */
/* ------------------------------------------------------------------ */

#define ERND MPFR_RNDN

typedef struct { mpfr_t re, im; } ecx;

static void ecx_init(ecx* z, mpfr_prec_t p) { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void ecx_clear(ecx* z)               { mpfr_clear(z->re);    mpfr_clear(z->im);    }
static void ecx_set(ecx* d, const ecx* s)   { mpfr_set(d->re, s->re, ERND); mpfr_set(d->im, s->im, ERND); }

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

/* out = exp(z) = e^re (cos im + i sin im). */
static void ecx_exp(ecx* out, const ecx* z, mpfr_prec_t p) {
    mpfr_t ea, c, s;
    mpfr_inits2(p, ea, c, s, (mpfr_ptr)0);
    mpfr_exp(ea, z->re, ERND);
    mpfr_sin_cos(s, c, z->im, ERND);
    mpfr_mul(out->re, ea, c, ERND);
    mpfr_mul(out->im, ea, s, ERND);
    mpfr_clears(ea, c, s, (mpfr_ptr)0);
}

/* erf(z) for complex z into already-init'd `out` (precision wp).
 * Returns false if the series fails to converge within the cap. */
static bool ecx_erf_series(ecx* out, const ecx* z, mpfr_prec_t wp, double z2abs) {
    ecx z2, two_z2, t, s, tmp;
    ecx_init(&z2, wp); ecx_init(&two_z2, wp);
    ecx_init(&t, wp); ecx_init(&s, wp); ecx_init(&tmp, wp);

    ecx_mul(&z2, z, z, wp);                  /* z^2 */
    mpfr_mul_2ui(two_z2.re, z2.re, 1, ERND); /* 2 z^2 */
    mpfr_mul_2ui(two_z2.im, z2.im, 1, ERND);
    ecx_set(&t, z);                          /* t_0 = z */
    ecx_set(&s, z);                          /* S    = z */

    mpfr_t mag, smag, eps;
    mpfr_inits2(wp, mag, smag, eps, (mpfr_ptr)0);
    mpfr_set_ui(eps, 1, ERND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 8 ? wp - 4 : 1), ERND); /* 2^-(wp-4) */

    unsigned long cap = (unsigned long)(2.0 * z2abs) + (unsigned long)wp + 1000;
    bool ok = false;
    for (unsigned long n = 1; n <= cap; n++) {
        ecx_mul(&tmp, &t, &two_z2, wp);
        ecx_div_d(&t, &tmp, (double)(2 * n + 1));   /* t_n = t_{n-1} (2 z^2)/(2n+1) */
        ecx_add(&s, &s, &t);
        if ((double)n > z2abs) {                    /* past the peak */
            ecx_abs(mag, &t);
            ecx_abs(smag, &s);
            mpfr_mul(smag, smag, eps, ERND);
            if (mpfr_cmp(mag, smag) < 0) { ok = true; break; }
        }
    }

    if (ok) {
        /* out = (2/sqrt(pi)) e^-z^2 S. */
        ecx nz2, ez;
        ecx_init(&nz2, wp); ecx_init(&ez, wp);
        mpfr_neg(nz2.re, z2.re, ERND);
        mpfr_neg(nz2.im, z2.im, ERND);
        ecx_exp(&ez, &nz2, wp);                     /* e^-z^2 */
        ecx_mul(out, &ez, &s, wp);

        mpfr_t coeff;
        mpfr_init2(coeff, wp);
        mpfr_const_pi(coeff, ERND);
        mpfr_sqrt(coeff, coeff, ERND);
        mpfr_ui_div(coeff, 2, coeff, ERND);         /* 2/sqrt(pi) */
        mpfr_mul(out->re, out->re, coeff, ERND);
        mpfr_mul(out->im, out->im, coeff, ERND);
        mpfr_clear(coeff);
        ecx_clear(&nz2); ecx_clear(&ez);
    }

    mpfr_clears(mag, smag, eps, (mpfr_ptr)0);
    ecx_clear(&z2); ecx_clear(&two_z2);
    ecx_clear(&t); ecx_clear(&s); ecx_clear(&tmp);
    return ok;
}

/* erfi(z) = -i erf(i z) for complex z into already-init'd `out`. */
static bool ecx_erfi(ecx* out, const ecx* z, mpfr_prec_t wp, double z2abs) {
    ecx iz, e;
    ecx_init(&iz, wp); ecx_init(&e, wp);
    mpfr_neg(iz.re, z->im, ERND);   /* Re(i z) = -Im(z) */
    mpfr_set(iz.im, z->re, ERND);   /* Im(i z) =  Re(z) */
    /* |i z|^2 = |z|^2, so z2abs carries over unchanged. */
    bool ok = ecx_erf_series(&e, &iz, wp, z2abs);
    if (ok) {
        /* -i (e.re + i e.im) = e.im - i e.re. */
        mpfr_set(out->re, e.im, ERND);
        mpfr_neg(out->im, e.re, ERND);
    }
    ecx_clear(&iz); ecx_clear(&e);
    return ok;
}

/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool erfi_set_mpfr(mpfr_t out, const Expr* e) {
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

/* Build a complex result from (re, im) at out_prec: machine precision (<= 53)
 * yields Real parts, higher yields MPFR. make_complex drops a zero imaginary
 * part to a bare real. */
static Expr* erfi_complex_result(const mpfr_t re, const mpfr_t im, mpfr_prec_t out_prec) {
    Expr *rr, *ii;
    if (out_prec <= 53) {
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

/* Build a real result at out_prec: Real if out_prec <= 53, else MPFR. */
static Expr* erfi_real_result(const mpfr_t v, mpfr_prec_t out_prec) {
    if (out_prec <= 53) return expr_new_real(mpfr_get_d(v, ERND));
    Expr* out = expr_new_mpfr_bits(out_prec);
    mpfr_set(out->data.mpfr, v, ERND);
    return out;
}

/* erfi(x) for a real MPFR x via the all-positive Maclaurin series
 *   erfi(x) = (2/sqrt(pi)) Sum_{n>=0} x^(2n+1)/(n!(2n+1)),
 *   u_0 = x,  u_n = u_{n-1} x^2 (2n-1)/(n(2n+1)).
 * For real x every term shares x's sign, so the partial sums climb
 * monotonically to the result -- no cancellation. Returns an Expr at out_prec,
 * or NULL on non-convergence. */
static Expr* erfi_mpfr_real(const mpfr_t x, mpfr_prec_t out_prec) {
    if (mpfr_zero_p(x)) {
        mpfr_t z; mpfr_init2(z, out_prec <= 53 ? 53 : out_prec);
        mpfr_set_zero(z, mpfr_signbit(x) ? -1 : 1);
        Expr* out = erfi_real_result(z, out_prec);
        mpfr_clear(z);
        return out;
    }

    double xd = mpfr_get_d(x, ERND);
    double x2abs = xd * xd;
    /* All-positive sum: no catastrophic cancellation, only accumulated
     * rounding over the iteration count -- a flat 64-bit guard covers it. */
    mpfr_prec_t wp = out_prec + 64;

    mpfr_t x2, u, s, mag, smag, eps;
    mpfr_inits2(wp, x2, u, s, mag, smag, eps, (mpfr_ptr)0);
    mpfr_set(u, x, ERND);                /* u_0 = x */
    mpfr_set(s, x, ERND);                /* S    = x */
    mpfr_mul(x2, u, u, ERND);            /* x^2 */
    mpfr_set_ui(eps, 1, ERND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 8 ? wp - 4 : 1), ERND); /* 2^-(wp-4) */

    unsigned long cap = (unsigned long)(2.0 * x2abs) + (unsigned long)wp + 1000;
    bool ok = false;
    for (unsigned long n = 1; n <= cap; n++) {
        mpfr_mul(u, u, x2, ERND);
        mpfr_mul_ui(u, u, 2 * n - 1, ERND);
        mpfr_div_ui(u, u, n, ERND);
        mpfr_div_ui(u, u, 2 * n + 1, ERND);   /* u_n = u_{n-1} x^2 (2n-1)/(n(2n+1)) */
        mpfr_add(s, s, u, ERND);
        if ((double)n > x2abs) {              /* past the term peak */
            mpfr_abs(mag, u, ERND);
            mpfr_abs(smag, s, ERND);
            mpfr_mul(smag, smag, eps, ERND);
            if (mpfr_cmp(mag, smag) < 0) { ok = true; break; }
        }
    }

    Expr* out = NULL;
    if (ok) {
        mpfr_t coeff;
        mpfr_init2(coeff, wp);
        mpfr_const_pi(coeff, ERND);
        mpfr_sqrt(coeff, coeff, ERND);
        mpfr_ui_div(coeff, 2, coeff, ERND);   /* 2/sqrt(pi) */
        mpfr_mul(s, s, coeff, ERND);
        mpfr_clear(coeff);
        if (!mpfr_nan_p(s) && !mpfr_inf_p(s))
            out = erfi_real_result(s, out_prec);
    }

    mpfr_clears(x2, u, s, mag, smag, eps, (mpfr_ptr)0);
    return out;
}

/* Evaluate erfi(z) for a numeric complex `arg` (Complex[..] with at least one
 * inexact part) at out_prec bits, via the MPFR series. Returns the result, or
 * NULL if the parts are not numeric / the series diverges. */
static Expr* erfi_mpfr_complex(Expr* re, Expr* im, mpfr_prec_t out_prec) {
    double red = 0.0, imd = 0.0;
    (void)erfi_to_double(re, &red);
    (void)erfi_to_double(im, &imd);
    double z2abs = red * red + imd * imd;
    long guard = 64 + (long)(z2abs / M_LN2);
    mpfr_prec_t wp = out_prec + (mpfr_prec_t)guard;

    ecx z, g;
    ecx_init(&z, wp); ecx_init(&g, wp);
    Expr* out = NULL;
    if (erfi_set_mpfr(z.re, re) && erfi_set_mpfr(z.im, im) &&
        ecx_erfi(&g, &z, wp, z2abs) &&
        !mpfr_nan_p(g.re) && !mpfr_nan_p(g.im) &&
        !mpfr_inf_p(g.re) && !mpfr_inf_p(g.im)) {
        out = erfi_complex_result(g.re, g.im, out_prec);
    }
    ecx_clear(&z); ecx_clear(&g);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Erfi[z]                                                            */
/* ------------------------------------------------------------------ */

static Expr* erfi_one_arg(Expr* arg) {
    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return expr_new_integer(0);                       /* Erfi[0] = 0 */
    if (erfi_is_symbol(arg, "Infinity"))      return expr_new_symbol(SYM_Infinity);
    if (erfi_is_neg_infinity(arg))
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2));
    if (erfi_is_symbol(arg, "ComplexInfinity")) return expr_new_symbol(SYM_ComplexInfinity);
    if (erfi_is_symbol(arg, "Indeterminate"))   return expr_new_symbol(SYM_Indeterminate);
    {
        int s = erfi_directed_imag_infinity(arg);
        if (s != 0)                                       /* +-I Infinity -> +-I */
            return make_complex(expr_new_integer(0), expr_new_integer(s));
    }

    /* 2. Machine real. */
    if (arg->type == EXPR_REAL) {
#ifdef USE_MPFR
        mpfr_t x;
        mpfr_init2(x, 53);
        mpfr_set_d(x, arg->data.real, ERND);
        Expr* out = erfi_mpfr_real(x, 53);
        mpfr_clear(x);
        if (out) return out;
#else
        double v;
        if (erfi_machine_real(arg->data.real, &v)) return expr_new_real(v);
#endif
    }

#ifdef USE_MPFR
    /* 3. Arbitrary-precision real. */
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        Expr* out = erfi_mpfr_real(arg->data.mpfr, prec);
        if (out) return out;
    }
#endif

    /* 4. Complex argument (Complex[..] with an inexact part). */
    {
        Expr *re, *im;
        if (is_complex(arg, &re, &im) && (erfi_is_inexact(re) || erfi_is_inexact(im))) {
#ifdef USE_MPFR
            long bits = numeric_min_inexact_bits(arg);
            mpfr_prec_t out_prec = (bits && bits > 53) ? (mpfr_prec_t)bits : 53;
            Expr* out = erfi_mpfr_complex(re, im, out_prec);
            if (out) return out;
#else
            double rr, ii;
            if (erfi_to_double(re, &rr) && erfi_to_double(im, &ii)) {
                double complex v;
                if (erfi_machine_complex(rr + ii * I, &v)) {
                    if (cimag(v) == 0.0) return expr_new_real(creal(v));
                    return make_complex(expr_new_real(creal(v)),
                                        expr_new_real(cimag(v)));
                }
            }
#endif
        }
    }

    /* 5. Odd symmetry: Erfi[-x] = -Erfi[x] for a symbolic negative-leading
     *    argument (numeric negatives are already handled above). */
    if (erfi_is_neg_leading_times(arg)) {
        Expr* pos = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1), expr_copy(arg) }, 2));
        Expr* inner = expr_new_function(expr_new_symbol(SYM_Erfi), &pos, 1);
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1), inner }, 2));
    }

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

/* Mathematica-compatible argx diagnostic; returns NULL so the evaluator
 * leaves the call unevaluated. */
static Expr* erfi_emit_argx(size_t argc) {
    fprintf(stderr,
            "Erfi::argx: Erfi called with %zu argument%s; "
            "1 argument is expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_erfi(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return erfi_one_arg(args[0]);
    return erfi_emit_argx(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void erfi_init(void) {
    symtab_add_builtin("Erfi", builtin_erfi);
    symtab_get_def("Erfi")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

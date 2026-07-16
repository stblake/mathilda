/* Mathilda -- the hyperbolic cosine integral
 *   Chi(z) = EulerGamma + Log[z] + Int_0^z (Cosh[t] - 1)/t dt.
 *
 *   CoshIntegral[z]   Chi(z)
 *
 * Chi has a logarithmic singularity at z = 0 and a branch cut running along the
 * negative real axis (-Infinity, 0].  It is the imaginary-axis sibling of the
 * cosine integral: Chi(z) = Ci(i z) - i Pi/2.  Evaluation is layered so each
 * argument takes the cheapest route:
 *
 *   exact special values     ->  -Infinity, Infinity, +-I Pi/2, Indeterminate
 *   machine real             ->  MPFR series/asymptotic at 53 bits
 *   arbitrary real           ->  the same, at the input precision
 *   complex (any precision)  ->  the ncpx series/asymptotic with guard bits
 *   everything else          ->  stays symbolic (return NULL)
 *
 * The convergent Maclaurin series (the trig series with the alternating sign
 * removed) is
 *
 *   Chi(z) = EulerGamma + Log(z) + Sum_{k>=1} z^(2k) / (2k (2k)!),
 *
 * valid on the principal branch for all z != 0 -- the principal Log(z) supplies
 * the correct +-i Pi jump across the cut, so no folding is needed for the
 * convergent path.  Every series term is positive, so on the real axis there is
 * NO catastrophic cancellation and only a small fixed guard is needed.  For large
 * |z| the convergent series is infeasible; there we use the asymptotic expansion
 *
 *   Chi(z) ~ sinh(z) F(z) + cosh(z) G(z)   (+ Stokes constant, see below),
 *     F(z) ~ Sum (2k)!   / z^(2k+1) = 1/z + 2!/z^3 + ...,
 *     G(z) ~ Sum (2k+1)! / z^(2k+2) = 1/z^2 + 3!/z^4 + ...,
 *
 * (the same F, G SinhIntegral computes; only the combination differs -- Shi uses
 * cosh F + sinh G).  The bare part sinh F + cosh G is even and real on the
 * positive real axis (no constant).  Off the real axis, matching Mathematica's
 * Series[CoshIntegral[z],{z,Infinity,k}], a piecewise Stokes term restores the
 * principal branch:
 *
 *   Chi(z) = B2(z) + i K,   B2(z) = sinh(z) F(z) + cosh(z) G(z),
 *     K = -Pi/2                       for Im(z) < 0,
 *     K = Pi sgn+(Re z) - Pi/2        for Im(z) > 0   (sgn+(0) = +1, from above).
 *
 * (On the imaginary axis Im(z) > 0, Re(z) = 0 gives K = +Pi/2, so
 * Chi(+- i Infinity) = +- i Pi/2.)  A negative real x is handled by the real
 * path as the from-above branch value Complex[Chi(|x|), Pi], matching the log
 * jump; and Chi(-Infinity) -> Infinity (the real part dominates).
 *
 * Machine-real results that overflow a C double (e.g. CoshIntegral[10.^6] ~
 * 1.5*10^434288) are emitted as a 53-bit MPFR real, whose exponent range easily
 * holds them, rather than as an inf-valued double.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "coshintegral.h"
#include "sym_names.h"

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "arithmetic.h"        /* is_rational, make_complex, is_complex */
#include "numeric.h"           /* numeric_min_inexact_bits */
#include "attr.h"
#include "eval.h"              /* eval_and_free */
#include "expr.h"
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#include "numeric_complex.h"   /* ncpx toolkit */
#define SRND MPFR_RNDN
#endif

/* M_LN2 is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

/* Pi and Pi/2 to double precision, for the non-MPFR fallback paths. */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

/* EulerGamma to double precision, for the non-MPFR fallback paths. */
#define CHI_EULER_GAMMA 0.57721566490153286061

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

static bool chi_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

static bool chi_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

static bool chi_to_double(const Expr* e, double* out) {
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
static bool chi_is_neg_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!chi_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == -1 && chi_is_symbol(b, "Infinity"))
        return true;
    if (b->type == EXPR_INTEGER && b->data.integer == -1 && chi_is_symbol(a, "Infinity"))
        return true;
    return false;
}

/* Match Times[Complex[0, s], Infinity] (s = +-1): +-I Infinity. Returns
 * +1 for I Infinity, -1 for -I Infinity, 0 if not of this form. */
static int chi_directed_imag_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return 0;
    if (!chi_is_symbol(e->data.function.head, "Times")) return 0;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    Expr *cplx = NULL, *inf = NULL;
    if (chi_is_symbol(b, "Infinity")) { cplx = a; inf = b; }
    else if (chi_is_symbol(a, "Infinity")) { cplx = b; inf = a; }
    if (!inf) return 0;
    Expr *re, *im;
    if (!is_complex(cplx, &re, &im)) return 0;
    if (!(re->type == EXPR_INTEGER && re->data.integer == 0)) return 0;
    if (im->type == EXPR_INTEGER && im->data.integer == 1)  return 1;
    if (im->type == EXPR_INTEGER && im->data.integer == -1) return -1;
    return 0;
}

/* Build -Infinity = Times[-1, Infinity]. */
static Expr* chi_make_neg_infinity(void) {
    return expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2);
}

/* Build +-I Pi/2 = Times[Complex[0, s], Rational[1, 2], Pi] for s = +-1. */
static Expr* chi_make_directed_imag_half_pi(int s) {
    Expr* imag = make_complex(expr_new_integer(0), expr_new_integer(s));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ imag, make_rational(1, 2), expr_new_symbol(SYM_Pi) }, 3));
}

#ifndef USE_MPFR
/* ------------------------------------------------------------------ */
/* Machine-precision fallback paths (USE_MPFR=0 builds)               */
/* ------------------------------------------------------------------ */

/* Chi(a) for a machine real a > 0 via convergent series / asymptotic. */
static bool chi_machine_pos(double a, double* out) {
    double v;
    if (a <= 40.0) {
        double a2 = a * a, term = a2 / 4.0, sum = a2 / 4.0;     /* k = 1 term */
        for (int k = 2; k <= 100000; k++) {
            term *= a2 * (double)(k - 1)
                    / ((double)k * (2 * k) * (2 * k - 1));
            sum += term;
            if ((double)(2 * k) > a && fabs(term) <= 1e-18 * (fabs(sum) + 1.0)) break;
        }
        v = CHI_EULER_GAMMA + log(a) + sum;
    } else {
        double f = 1.0 / a, g = 1.0 / (a * a), a2 = a * a;
        double tf = f, tg = g, pf = f, pg = g;
        for (int k = 1; k <= (int)a; k++) {
            tf *= (double)(2 * k) * (2 * k - 1) / a2;
            if (fabs(tf) >= fabs(pf)) break;
            f += tf; pf = tf;
        }
        for (int k = 1; k <= (int)a; k++) {
            tg *= (double)(2 * k + 1) * (2 * k) / a2;
            if (fabs(tg) >= fabs(pg)) break;
            g += tg; pg = tg;
        }
        v = sinh(a) * f + cosh(a) * g;
    }
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

/* Chi(x) for a machine real x != 0. For x < 0 the from-above branch value is
 * Chi(|x|) + i Pi, so (*im) carries Pi. */
static bool chi_machine_real(double x, double* re, double* im) {
    double v;
    if (!chi_machine_pos(fabs(x), &v)) return false;
    *re = v;
    *im = (x < 0.0) ? M_PI : 0.0;
    return true;
}

/* Chi(z) for a machine double-complex z via the convergent series. */
static bool chi_machine_complex(double complex z, double complex* out) {
    double complex z2 = z * z, term = z2 / 4.0, sum = z2 / 4.0;   /* k = 1 term */
    double zabs = cabs(z);
    for (int k = 2; k <= 200000; k++) {
        term *= z2 * (double)(k - 1) / ((double)k * (2 * k) * (2 * k - 1));
        sum += term;
        if ((double)(2 * k) > zabs && cabs(term) <= 1e-17 * (cabs(sum) + 1.0)) break;
    }
    double complex v = CHI_EULER_GAMMA + clog(z) + sum;   /* clog: principal */
    if (!isfinite(creal(v)) || !isfinite(cimag(v))) return false;
    *out = v;
    return true;
}
#endif /* !USE_MPFR */

#ifdef USE_MPFR
/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool chi_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, SRND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          SRND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        SRND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          SRND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, SRND);
        mpfr_div_si(out, out, (long)d, SRND);
        return true;
    }
    return false;
}

/* Emit one real component `v` at out_prec. For out_prec <= 53 the value normally
 * becomes a machine Real; but when it overflows the double range (Chi grows like
 * e^|x|), fall back to a 53-bit MPFR real so the huge exponent is preserved
 * instead of collapsing to inf. */
static Expr* chi_emit_component(const mpfr_t v, mpfr_prec_t out_prec) {
    if (out_prec <= 53) {
        double d = mpfr_get_d(v, SRND);
        if (isfinite(d) || mpfr_inf_p(v) || mpfr_zero_p(v))
            return expr_new_real(d);
        Expr* e = expr_new_mpfr_bits(53);
        mpfr_set(e->data.mpfr, v, SRND);
        return e;
    }
    Expr* e = expr_new_mpfr_bits(out_prec);
    mpfr_set(e->data.mpfr, v, SRND);
    return e;
}

/* -------------------- real paths -------------------- */

/* Convergent even series Sum_{k>=1} a^(2k)/(2k (2k)!) into `sum` (precision wp),
 * for a > 0. All terms positive. Stops once past 2k > a and the term has decayed
 * below 2^-(wp-8) of the largest term seen. */
static void chi_real_series(mpfr_t sum, const mpfr_t a, mpfr_prec_t wp, double aabs) {
    mpfr_t term, a2, mag, peak, thr;
    mpfr_inits2(wp, term, a2, mag, peak, thr, (mpfr_ptr)0);
    mpfr_mul(a2, a, a, SRND);
    mpfr_div_ui(term, a2, 4, SRND);       /* k = 1 term = a^2/4 */
    mpfr_set(sum, term, SRND);
    mpfr_abs(peak, term, SRND);
    long shift = (long)(wp > 8 ? wp - 8 : 1);
    unsigned long cap = (unsigned long)(2.0 * aabs) + (unsigned long)wp + 1000;
    for (unsigned long k = 2; k <= cap; k++) {
        /* term *= a^2 (k-1) / (k (2k)(2k-1)). */
        mpfr_mul(term, term, a2, SRND);
        mpfr_mul_ui(term, term, k - 1, SRND);
        mpfr_div_ui(term, term, k, SRND);
        mpfr_div_ui(term, term, 2 * k, SRND);
        mpfr_div_ui(term, term, 2 * k - 1, SRND);
        mpfr_add(sum, sum, term, SRND);
        mpfr_abs(mag, term, SRND);
        if (mpfr_cmp(mag, peak) > 0) mpfr_set(peak, mag, SRND);
        if ((double)(2 * k) > aabs) {
            mpfr_mul_2si(thr, peak, -shift, SRND);
            if (mpfr_cmp(mag, thr) < 0) break;
        }
    }
    mpfr_clears(term, a2, mag, peak, thr, (mpfr_ptr)0);
}

/* Asymptotic Chi(a) = sinh(a) F(a) + cosh(a) G(a) into `out` (precision wp), for
 * a > 0. Each of F, G (all-positive) is summed to its smallest term. */
static void chi_real_asymptotic(mpfr_t out, const mpfr_t a, mpfr_prec_t wp) {
    mpfr_t f, g, tf, tg, a2, mag, prev, thr, ca, sa;
    mpfr_inits2(wp, f, g, tf, tg, a2, mag, prev, thr, ca, sa, (mpfr_ptr)0);
    mpfr_mul(a2, a, a, SRND);
    mpfr_set_ui(thr, 1, SRND);
    mpfr_mul_2si(thr, thr, -(long)wp, SRND);          /* 2^-wp */

    /* F = Sum (2k)!/a^(2k+1), t_0 = 1/a. */
    mpfr_ui_div(tf, 1, a, SRND);
    mpfr_set(f, tf, SRND);
    mpfr_abs(prev, tf, SRND);
    for (unsigned long k = 1; k <= 100000; k++) {
        mpfr_mul_si(tf, tf, (long)((2 * k) * (2 * k - 1)), SRND);
        mpfr_div(tf, tf, a2, SRND);
        mpfr_abs(mag, tf, SRND);
        if (mpfr_cmp(mag, prev) >= 0) break;          /* terms growing */
        mpfr_add(f, f, tf, SRND);
        mpfr_set(prev, mag, SRND);
        if (mpfr_cmp(mag, thr) < 0) break;
    }
    /* G = Sum (2k+1)!/a^(2k+2), t_0 = 1/a^2. */
    mpfr_ui_div(tg, 1, a2, SRND);
    mpfr_set(g, tg, SRND);
    mpfr_abs(prev, tg, SRND);
    for (unsigned long k = 1; k <= 100000; k++) {
        mpfr_mul_si(tg, tg, (long)((2 * k + 1) * (2 * k)), SRND);
        mpfr_div(tg, tg, a2, SRND);
        mpfr_abs(mag, tg, SRND);
        if (mpfr_cmp(mag, prev) >= 0) break;
        mpfr_add(g, g, tg, SRND);
        mpfr_set(prev, mag, SRND);
        if (mpfr_cmp(mag, thr) < 0) break;
    }

    mpfr_sinh_cosh(sa, ca, a, SRND);
    mpfr_mul(sa, sa, f, SRND);                         /* sinh(a) F */
    mpfr_mul(ca, ca, g, SRND);                         /* cosh(a) G */
    mpfr_add(out, sa, ca, SRND);                       /* sinh F + cosh G */
    mpfr_clears(f, g, tf, tg, a2, mag, prev, thr, ca, sa, (mpfr_ptr)0);
}

/* Chi(a) for a > 0 into already-init'd `val` (precision wp = mpfr_get_prec(val)). */
static void chi_real_core(mpfr_t val, const mpfr_t a, double aabs, mpfr_prec_t out_prec) {
    if (aabs > (double)(out_prec + 24) * M_LN2 + 8.0) {
        chi_real_asymptotic(val, a, mpfr_get_prec(val));
    } else {
        mpfr_t euler, lna, s;
        mpfr_prec_t wp = mpfr_get_prec(val);
        mpfr_inits2(wp, euler, lna, s, (mpfr_ptr)0);
        chi_real_series(s, a, wp, aabs);
        mpfr_const_euler(euler, SRND);
        mpfr_log(lna, a, SRND);
        mpfr_add(val, euler, lna, SRND);
        mpfr_add(val, val, s, SRND);                   /* gamma + ln a + series */
        mpfr_clears(euler, lna, s, (mpfr_ptr)0);
    }
}

/* Build a real-or-complex result (re, im) at out_prec, overflow-aware per
 * component. make_complex drops a zero imaginary part to a bare real. */
static Expr* chi_complex_result(const mpfr_t re, const mpfr_t im, mpfr_prec_t out_prec) {
    return make_complex(chi_emit_component(re, out_prec),
                        chi_emit_component(im, out_prec));
}

/* Chi(x) for a numeric real `arg` at out_prec bits. x > 0 yields a bare real;
 * x < 0 yields the from-above branch value Complex[Chi(|x|), Pi]; x == 0 yields
 * -Infinity. */
static Expr* chi_mpfr_real(const Expr* arg, mpfr_prec_t out_prec) {
    mpfr_t x;
    mpfr_init2(x, out_prec + 16);
    if (!chi_set_mpfr(x, arg)) { mpfr_clear(x); return NULL; }
    int sgn = mpfr_sgn(x);
    if (sgn == 0) { mpfr_clear(x); return chi_make_neg_infinity(); }

    double aabs = fabs(mpfr_get_d(x, SRND));
    /* All-positive series -> no cancellation guard; fixed 64-bit guard. */
    mpfr_prec_t wp = out_prec + 64;

    mpfr_t a, val;
    mpfr_inits2(wp, a, val, (mpfr_ptr)0);
    mpfr_abs(a, x, SRND);
    chi_real_core(val, a, aabs, out_prec);

    Expr* out;
    if (sgn > 0) {
        out = chi_emit_component(val, out_prec);
    } else {
        mpfr_t pi;
        mpfr_init2(pi, wp);
        mpfr_const_pi(pi, SRND);                        /* from-above jump: +Pi */
        out = chi_complex_result(val, pi, out_prec);
        mpfr_clear(pi);
    }
    mpfr_clears(a, val, (mpfr_ptr)0);
    mpfr_clear(x);
    return out;
}

/* -------------------- complex paths (ncpx) -------------------- */

/* Multiply the components of `z` in place by the integer m. */
static void chi_ncpx_mul_si(ncpx* z, long m) {
    mpfr_mul_si(z->re, z->re, m, SRND);
    mpfr_mul_si(z->im, z->im, m, SRND);
}

/* Divide the components of `z` in place by the positive integer m. */
static void chi_ncpx_div_ui(ncpx* z, unsigned long m) {
    mpfr_div_ui(z->re, z->re, m, SRND);
    mpfr_div_ui(z->im, z->im, m, SRND);
}

/* cosh(w) and sinh(w) into cw, sw via the exponential (ncpx has no hyperbolic
 * primitives): cosh = (e^w + e^-w)/2, sinh = (e^w - e^-w)/2. */
static void chi_ncpx_cosh_sinh(ncpx* cw, ncpx* sw, const ncpx* w, mpfr_prec_t wp) {
    ncpx ew, enw, nw;
    ncpx_init(&ew, wp); ncpx_init(&enw, wp); ncpx_init(&nw, wp);
    ncpx_exp(&ew, w, wp);
    ncpx_neg(&nw, w);
    ncpx_exp(&enw, &nw, wp);
    ncpx_add(cw, &ew, &enw);
    ncpx_sub(sw, &ew, &enw);
    mpfr_div_2ui(cw->re, cw->re, 1, SRND);
    mpfr_div_2ui(cw->im, cw->im, 1, SRND);
    mpfr_div_2ui(sw->re, sw->re, 1, SRND);
    mpfr_div_2ui(sw->im, sw->im, 1, SRND);
    ncpx_clear(&ew); ncpx_clear(&enw); ncpx_clear(&nw);
}

/* Convergent series Chi(z) = gamma + Log(z) + Sum z^(2k)/(2k (2k)!) into `out`,
 * valid on the principal branch for any z != 0. All terms positive. Returns
 * false if the even series fails to converge within the cap. */
static bool chi_ncpx_series(ncpx* out, const ncpx* z, mpfr_prec_t wp, double zabs) {
    ncpx term, z2, sum, lg;
    ncpx_init(&term, wp); ncpx_init(&z2, wp);
    ncpx_init(&sum, wp);  ncpx_init(&lg, wp);
    ncpx_mul(&z2, z, z, wp);
    ncpx_set(&term, &z2);                  /* k = 1 term = z^2/4 */
    chi_ncpx_div_ui(&term, 4);
    ncpx_set(&sum, &term);

    mpfr_t mag, peak, thr, euler;
    mpfr_inits2(wp, mag, peak, thr, euler, (mpfr_ptr)0);
    ncpx_abs(peak, &term);
    long shift = (long)(wp > 8 ? wp - 8 : 1);
    unsigned long cap = (unsigned long)(2.0 * zabs) + (unsigned long)wp + 1000;
    bool ok = false;
    for (unsigned long k = 2; k <= cap; k++) {
        /* term *= z^2 (k-1) / (k (2k)(2k-1)). */
        ncpx_mul(&term, &term, &z2, wp);
        chi_ncpx_mul_si(&term, (long)(k - 1));
        chi_ncpx_div_ui(&term, k);
        chi_ncpx_div_ui(&term, 2 * k);
        chi_ncpx_div_ui(&term, 2 * k - 1);
        ncpx_add(&sum, &sum, &term);
        ncpx_abs(mag, &term);
        if (mpfr_cmp(mag, peak) > 0) mpfr_set(peak, mag, SRND);
        if ((double)(2 * k) > zabs) {
            mpfr_mul_2si(thr, peak, -shift, SRND);
            if (mpfr_cmp(mag, thr) < 0) { ok = true; break; }
        }
    }
    if (ok) {
        ncpx_log(&lg, z, wp);                           /* principal Log(z) */
        mpfr_const_euler(euler, SRND);
        mpfr_add(out->re, euler, lg.re, SRND);
        mpfr_set(out->im, lg.im, SRND);
        ncpx_add(out, out, &sum);                       /* gamma + Log(z) + series */
    }
    mpfr_clears(mag, peak, thr, euler, (mpfr_ptr)0);
    ncpx_clear(&term); ncpx_clear(&z2); ncpx_clear(&sum); ncpx_clear(&lg);
    return ok;
}

/* Sum one auxiliary asymptotic series into `acc` (already = t_0), advancing the
 * running term `t` by the factor (a)(b)/w2 each step to its smallest term.
 * `odd` selects the F-series (base 2k, 2k-1) vs G-series (2k+1, 2k). All terms
 * positive (no sign alternation). */
static void chi_ncpx_aux(ncpx* acc, ncpx* t, const ncpx* w2, mpfr_prec_t wp,
                         mpfr_t thr, int odd) {
    mpfr_t mag, prev;
    mpfr_inits2(wp, mag, prev, (mpfr_ptr)0);
    ncpx_abs(prev, t);
    for (unsigned long k = 1; k <= 100000; k++) {
        long a = odd ? (long)(2 * k) : (long)(2 * k + 1);
        long b = odd ? (long)(2 * k - 1) : (long)(2 * k);
        chi_ncpx_mul_si(t, a);
        chi_ncpx_mul_si(t, b);
        ncpx_div(t, t, w2, wp);
        ncpx_abs(mag, t);
        if (mpfr_cmp(mag, prev) >= 0) break;
        ncpx_add(acc, acc, t);
        mpfr_set(prev, mag, SRND);
        if (mpfr_cmp(mag, thr) < 0) break;
    }
    mpfr_clears(mag, prev, (mpfr_ptr)0);
}

/* Bare asymptotic B2(w) = sinh(w) F(w) + cosh(w) G(w) into `out`, for Re(w) >= 0.
 * B2 is even, so B2(z) = B2(w) after the fold; the caller adds the Stokes term. */
static void chi_ncpx_asymptotic(ncpx* out, const ncpx* w, mpfr_prec_t wp) {
    ncpx w2, one, f, g, tf, tg, cw, sw;
    ncpx_init(&w2, wp); ncpx_init(&one, wp);
    ncpx_init(&f, wp);  ncpx_init(&g, wp);
    ncpx_init(&tf, wp); ncpx_init(&tg, wp);
    ncpx_init(&cw, wp); ncpx_init(&sw, wp);
    mpfr_t thr;
    mpfr_init2(thr, wp);
    mpfr_set_ui(thr, 1, SRND);
    mpfr_mul_2si(thr, thr, -(long)wp, SRND);

    ncpx_mul(&w2, w, w, wp);
    ncpx_set_ui(&one, 1);
    ncpx_div(&tf, &one, w, wp);        /* 1/w    */
    ncpx_set(&f, &tf);
    chi_ncpx_aux(&f, &tf, &w2, wp, thr, 1);
    ncpx_div(&tg, &one, &w2, wp);      /* 1/w^2  */
    ncpx_set(&g, &tg);
    chi_ncpx_aux(&g, &tg, &w2, wp, thr, 0);

    chi_ncpx_cosh_sinh(&cw, &sw, w, wp);
    ncpx_mul(&sw, &sw, &f, wp);        /* sinh(w) F */
    ncpx_mul(&cw, &cw, &g, wp);        /* cosh(w) G */
    ncpx_add(out, &sw, &cw);           /* sinh F + cosh G */

    mpfr_clear(thr);
    ncpx_clear(&w2); ncpx_clear(&one); ncpx_clear(&f); ncpx_clear(&g);
    ncpx_clear(&tf); ncpx_clear(&tg); ncpx_clear(&cw); ncpx_clear(&sw);
}

/* Chi(z) for a numeric complex `arg` at out_prec bits. Moderate |z| uses the
 * principal convergent series (correct across the whole plane); large |z| uses
 * the right-half asymptotic with the piecewise Stokes constant. */
static Expr* chi_mpfr_complex(Expr* re, Expr* im, mpfr_prec_t out_prec) {
    double red = 0.0, imd = 0.0;
    (void)chi_to_double(re, &red);
    (void)chi_to_double(im, &imd);
    double zabs = sqrt(red * red + imd * imd);

    Expr* out = NULL;
    if (zabs > (double)(out_prec + 24) * M_LN2 + 8.0) {
        /* Asymptotic. Fold Re(z) < 0 to the right half (B2 is even, so the bare
         * part is unchanged); the Stokes constant K depends on z's location. */
        mpfr_prec_t wp = out_prec + 96;
        int fold = (red < 0.0);
        ncpx w, g;
        ncpx_init(&w, wp); ncpx_init(&g, wp);
        if (chi_set_mpfr(w.re, re) && chi_set_mpfr(w.im, im)) {
            if (fold) ncpx_neg(&w, &w);                 /* w = -z, Re(w) >= 0 */
            chi_ncpx_asymptotic(&g, &w, wp);            /* B2(z) = B2(w) */
            if (imd != 0.0) {
                mpfr_t pi, k;
                mpfr_init2(pi, wp);
                mpfr_init2(k, wp);
                mpfr_const_pi(pi, SRND);
                if (imd < 0.0) {
                    mpfr_div_2ui(k, pi, 1, SRND);       /* Pi/2 */
                    mpfr_neg(k, k, SRND);               /* K = -Pi/2 */
                } else {
                    /* K = Pi sgn+(Re z) - Pi/2, sgn+(0) = +1 (from above). */
                    mpfr_set(k, pi, SRND);
                    if (red < 0.0) mpfr_neg(k, k, SRND);
                    mpfr_t hp;
                    mpfr_init2(hp, wp);
                    mpfr_div_2ui(hp, pi, 1, SRND);
                    mpfr_sub(k, k, hp, SRND);
                    mpfr_clear(hp);
                }
                mpfr_add(g.im, g.im, k, SRND);          /* + i K */
                mpfr_clear(pi);
                mpfr_clear(k);
            }
            if (!mpfr_nan_p(g.re) && !mpfr_nan_p(g.im) &&
                !mpfr_inf_p(g.re) && !mpfr_inf_p(g.im))
                out = chi_complex_result(g.re, g.im, out_prec);
        }
        ncpx_clear(&w); ncpx_clear(&g);
        return out;
    }

    /* Convergent series: partial sums peak ~e^|z|; add |z|/ln2 guard bits to
     * absorb the cancellation exactly. */
    long guard = 64 + (long)(zabs / M_LN2);
    mpfr_prec_t wp = out_prec + (mpfr_prec_t)guard;
    ncpx z, g;
    ncpx_init(&z, wp); ncpx_init(&g, wp);
    if (chi_set_mpfr(z.re, re) && chi_set_mpfr(z.im, im)) {
        if (chi_ncpx_series(&g, &z, wp, zabs)) {
            if (!mpfr_nan_p(g.re) && !mpfr_nan_p(g.im) &&
                !mpfr_inf_p(g.re) && !mpfr_inf_p(g.im))
                out = chi_complex_result(g.re, g.im, out_prec);
        }
    }
    ncpx_clear(&z); ncpx_clear(&g);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* CoshIntegral[z]                                                    */
/* ------------------------------------------------------------------ */

static Expr* chi_one_arg(Expr* arg) {
    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return chi_make_neg_infinity();                      /* Chi[0] = -Inf */
    if (chi_is_symbol(arg, "Infinity"))         return expr_new_symbol(SYM_Infinity);
    if (chi_is_neg_infinity(arg))               return expr_new_symbol(SYM_Infinity);
    if (chi_is_symbol(arg, "ComplexInfinity"))  return expr_new_symbol(SYM_Indeterminate);
    if (chi_is_symbol(arg, "Indeterminate"))    return expr_new_symbol(SYM_Indeterminate);
    {
        int s = chi_directed_imag_infinity(arg);
        if (s != 0) return chi_make_directed_imag_half_pi(s);  /* +-I Pi/2 */
    }

    /* 2. Numeric real (machine and arbitrary precision). */
    if (arg->type == EXPR_REAL) {
#ifdef USE_MPFR
        return chi_mpfr_real(arg, 53);
#else
        double x = arg->data.real;
        if (x == 0.0) return chi_make_neg_infinity();
        double rr, ii;
        if (chi_machine_real(x, &rr, &ii)) {
            if (ii == 0.0) return expr_new_real(rr);
            return make_complex(expr_new_real(rr), expr_new_real(ii));
        }
        return NULL;
#endif
    }
#ifdef USE_MPFR
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        return chi_mpfr_real(arg, prec);
    }
#endif

    /* 3. Numeric complex (Complex[..] with an inexact part). */
    {
        Expr *re, *im;
        if (is_complex(arg, &re, &im) && (chi_is_inexact(re) || chi_is_inexact(im))) {
#ifdef USE_MPFR
            long bits = numeric_min_inexact_bits(arg);
            mpfr_prec_t out_prec = (bits && bits > 53) ? (mpfr_prec_t)bits : 53;
            Expr* out = chi_mpfr_complex(re, im, out_prec);
            if (out) return out;
#else
            double rr, ii;
            if (chi_to_double(re, &rr) && chi_to_double(im, &ii)) {
                double complex v;
                if (chi_machine_complex(rr + ii * I, &v)) {
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

static Expr* chi_emit_argx(size_t argc) {
    fprintf(stderr,
            "CoshIntegral::argx: CoshIntegral called with %zu arguments; "
            "1 argument is expected.\n",
            argc);
    return NULL;
}

Expr* builtin_coshintegral(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return chi_one_arg(args[0]);
    return chi_emit_argx(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void coshintegral_init(void) {
    symtab_add_builtin("CoshIntegral", builtin_coshintegral);
    symtab_get_def("CoshIntegral")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

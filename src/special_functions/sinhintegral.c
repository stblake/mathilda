/* Mathilda -- the hyperbolic sine integral  Shi(z) = Int_0^z Sinh[t]/t dt.
 *
 *   SinhIntegral[z]   Shi(z)
 *
 * Shi is entire and odd, with no branch cuts. It is the imaginary-axis sibling
 * of the sine integral: Shi(z) = -i Si(i z). Evaluation is layered so each kind
 * of argument takes the cheapest route:
 *
 *   exact special values     ->  0, +-Infinity, +-I Pi/2, Indeterminate
 *   machine real             ->  MPFR series/asymptotic at 53 bits
 *   arbitrary real           ->  the same, at the input precision
 *   complex (any precision)  ->  the ncpx series/asymptotic with guard bits
 *   everything else          ->  stays symbolic (return NULL)
 *
 * The convergent Maclaurin series (the trig series with the alternating sign
 * removed) is
 *
 *   Shi(z) = Sum_{k>=0} z^(2k+1) / ((2k+1) (2k+1)!),
 *
 * valid for all z. Unlike Si, every term is positive, so on the real axis there
 * is NO catastrophic cancellation -- the partial sums climb monotonically to the
 * O(e^|z|)-sized answer and only a small fixed guard is needed. (Complex inputs
 * still cancel like Si, so the complex paths keep the ~|z|/ln2 guard.) For large
 * |z| the convergent series is infeasible; there we use the asymptotic expansion
 *
 *   Shi(z) ~ cosh(z) F(z) + sinh(z) G(z)   (+ Stokes constant, see below),
 *     F(z) ~ Sum (2k)!   / z^(2k+1) = 1/z + 2!/z^3 + ...,
 *     G(z) ~ Sum (2k+1)! / z^(2k+2) = 1/z^2 + 3!/z^4 + ...,
 *
 * summed to the smallest term (optimal truncation). F, G are Si's f, g without
 * the (-1)^k. The odd-symmetry reduction Shi(-z) = -Shi(z) folds negative real /
 * left-half-plane inputs onto Re >= 0, keeping the asymptotic within its valid
 * sector. On the real axis the bare form above is exact (real, no constant); off
 * it, matching Mathematica's Series[SinhIntegral[z],{z,Infinity,k}], a Stokes
 * constant restores the analytic value:
 *
 *   Shi(z) = B(z) + i (Pi/2) sign(Im z),   B(z) = cosh(z) F(z) + sinh(z) G(z).
 *
 * (B is odd; on the imaginary axis B ~ O(1/z) vanishes and the i Pi/2 sign(Im z)
 * is all that survives, giving Shi(+- i Infinity) = +- i Pi/2.)
 *
 * Machine-real results that overflow a C double (e.g. SinhIntegral[10.^6] ~
 * 1.5*10^434288) are emitted as a 53-bit MPFR real, whose exponent range easily
 * holds them, rather than as an inf-valued double.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "sinhintegral.h"
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

/* Pi/2 to double precision, for the non-MPFR fallback paths. */
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

static bool shi_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

static bool shi_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

static bool shi_to_double(const Expr* e, double* out) {
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
static bool shi_is_neg_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!shi_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == -1 && shi_is_symbol(b, "Infinity"))
        return true;
    if (b->type == EXPR_INTEGER && b->data.integer == -1 && shi_is_symbol(a, "Infinity"))
        return true;
    return false;
}

/* Match Times[Complex[0, s], Infinity] (s = +-1): +-I Infinity. Returns
 * +1 for I Infinity, -1 for -I Infinity, 0 if not of this form. */
static int shi_directed_imag_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return 0;
    if (!shi_is_symbol(e->data.function.head, "Times")) return 0;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    Expr *cplx = NULL, *inf = NULL;
    if (shi_is_symbol(b, "Infinity")) { cplx = a; inf = b; }
    else if (shi_is_symbol(a, "Infinity")) { cplx = b; inf = a; }
    if (!inf) return 0;
    Expr *re, *im;
    if (!is_complex(cplx, &re, &im)) return 0;
    if (!(re->type == EXPR_INTEGER && re->data.integer == 0)) return 0;
    if (im->type == EXPR_INTEGER && im->data.integer == 1)  return 1;
    if (im->type == EXPR_INTEGER && im->data.integer == -1) return -1;
    return 0;
}

/* Build +-Infinity = Times[sign, Infinity] (sign = +-1). */
static Expr* shi_make_signed_infinity(int sign) {
    if (sign >= 0) return expr_new_symbol(SYM_Infinity);
    return expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2);
}

/* Build +-I Pi/2 = Times[Complex[0, s], Rational[1, 2], Pi] for s = +-1. */
static Expr* shi_make_directed_imag_half_pi(int s) {
    Expr* imag = make_complex(expr_new_integer(0), expr_new_integer(s));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ imag, make_rational(1, 2), expr_new_symbol(SYM_Pi) }, 3));
}

/* True if `arg` is Times[negative-literal, ...] -- the odd-symmetry trigger for
 * pulling a sign out of a symbolic argument. Mirrors sinintegral.c. */
static bool shi_is_neg_leading_times(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count < 2) return false;
    if (!shi_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    if (a->type == EXPR_INTEGER)  return a->data.integer < 0;
    int64_t n, d;
    if (is_rational(a, &n, &d))   return n < 0;
    if (a->type == EXPR_REAL)     return a->data.real < 0.0;
    return false;
}

#ifndef USE_MPFR
/* ------------------------------------------------------------------ */
/* Machine-precision fallback paths (USE_MPFR=0 builds)               */
/* ------------------------------------------------------------------ */

/* Shi(x) for a machine real x via the convergent series / asymptotic. Returns
 * false on overflow (|x| large: cosh/sinh exceed the double range). */
static bool shi_machine_real(double x, double* out) {
    double ax = fabs(x);
    double v;
    if (ax <= 40.0) {
        double term = ax, sum = ax, x2 = ax * ax;
        for (int k = 1; k <= 100000; k++) {
            term *= x2 * (2 * k - 1) / ((double)(2 * k) * (2 * k + 1) * (2 * k + 1));
            sum += term;
            if ((double)(2 * k) > ax && fabs(term) <= 1e-18 * (fabs(sum) + 1.0)) break;
        }
        v = sum;
    } else {
        double f = 1.0 / ax, g = 1.0 / (ax * ax), x2 = ax * ax;
        double tf = f, tg = g, pf = f, pg = g;
        for (int k = 1; k <= (int)ax; k++) {
            tf *= (double)(2 * k) * (2 * k - 1) / x2;
            if (fabs(tf) >= fabs(pf)) break;
            f += tf; pf = tf;
        }
        for (int k = 1; k <= (int)ax; k++) {
            tg *= (double)(2 * k + 1) * (2 * k) / x2;
            if (fabs(tg) >= fabs(pg)) break;
            g += tg; pg = tg;
        }
        v = cosh(ax) * f + sinh(ax) * g;
    }
    if (!isfinite(v)) return false;
    *out = (x < 0.0) ? -v : v;
    return true;
}

/* Shi(z) for a machine double-complex z via the convergent series. */
static bool shi_machine_complex(double complex z, double complex* out) {
    double complex w = (creal(z) < 0.0) ? -z : z;
    double zabs = cabs(w);
    double complex term = w, sum = w, w2 = w * w;
    for (int k = 1; k <= 200000; k++) {
        term *= w2 * (double)(2 * k - 1)
                / ((double)(2 * k) * (2 * k + 1) * (2 * k + 1));
        sum += term;
        if ((double)(2 * k) > zabs && cabs(term) <= 1e-17 * (cabs(sum) + 1.0)) break;
    }
    double complex v = (creal(z) < 0.0) ? -sum : sum;
    if (!isfinite(creal(v)) || !isfinite(cimag(v))) return false;
    *out = v;
    return true;
}
#endif /* !USE_MPFR */

#ifdef USE_MPFR
/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool shi_set_mpfr(mpfr_t out, const Expr* e) {
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
 * becomes a machine Real; but when it overflows the double range (Shi grows like
 * e^|x|), fall back to a 53-bit MPFR real so the huge exponent is preserved
 * instead of collapsing to inf. */
static Expr* shi_emit_component(const mpfr_t v, mpfr_prec_t out_prec) {
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

/* Convergent series Sum_{k>=0} x^(2k+1)/((2k+1)(2k+1)!) into `sum` (precision
 * wp), for x >= 0. All terms positive. Stops once past 2k > x and the term has
 * decayed below 2^-(wp-8) of the largest term seen. */
static void shi_real_series(mpfr_t sum, const mpfr_t x, mpfr_prec_t wp, double xabs) {
    mpfr_t term, x2, mag, peak, thr;
    mpfr_inits2(wp, term, x2, mag, peak, thr, (mpfr_ptr)0);
    mpfr_set(term, x, SRND);              /* k = 0 term = x        */
    mpfr_set(sum, x, SRND);
    mpfr_mul(x2, x, x, SRND);
    mpfr_abs(peak, x, SRND);
    long shift = (long)(wp > 8 ? wp - 8 : 1);
    unsigned long cap = (unsigned long)(2.0 * xabs) + (unsigned long)wp + 1000;
    for (unsigned long k = 1; k <= cap; k++) {
        mpfr_mul(term, term, x2, SRND);
        mpfr_mul_si(term, term, (long)(2 * k - 1), SRND);
        mpfr_div_si(term, term, (long)(2 * k), SRND);
        mpfr_div_si(term, term, (long)(2 * k + 1), SRND);
        mpfr_div_si(term, term, (long)(2 * k + 1), SRND);
        mpfr_add(sum, sum, term, SRND);
        mpfr_abs(mag, term, SRND);
        if (mpfr_cmp(mag, peak) > 0) mpfr_set(peak, mag, SRND);
        if ((double)(2 * k) > xabs) {
            mpfr_mul_2si(thr, peak, -shift, SRND);
            if (mpfr_cmp(mag, thr) < 0) break;
        }
    }
    mpfr_clears(term, x2, mag, peak, thr, (mpfr_ptr)0);
}

/* Asymptotic Shi(x) = cosh(x) F(x) + sinh(x) G(x) into `out` (precision wp), for
 * x > 0. Each of F, G (all-positive) is summed to its smallest term. */
static void shi_real_asymptotic(mpfr_t out, const mpfr_t x, mpfr_prec_t wp) {
    mpfr_t f, g, tf, tg, x2, mag, prev, thr, cx, sx;
    mpfr_inits2(wp, f, g, tf, tg, x2, mag, prev, thr, cx, sx, (mpfr_ptr)0);
    mpfr_mul(x2, x, x, SRND);
    mpfr_set_ui(thr, 1, SRND);
    mpfr_mul_2si(thr, thr, -(long)wp, SRND);          /* 2^-wp */

    /* F = Sum (2k)!/x^(2k+1), t_0 = 1/x. */
    mpfr_ui_div(tf, 1, x, SRND);
    mpfr_set(f, tf, SRND);
    mpfr_abs(prev, tf, SRND);
    for (unsigned long k = 1; k <= 100000; k++) {
        mpfr_mul_si(tf, tf, (long)((2 * k) * (2 * k - 1)), SRND);
        mpfr_div(tf, tf, x2, SRND);
        mpfr_abs(mag, tf, SRND);
        if (mpfr_cmp(mag, prev) >= 0) break;          /* terms growing */
        mpfr_add(f, f, tf, SRND);
        mpfr_set(prev, mag, SRND);
        if (mpfr_cmp(mag, thr) < 0) break;
    }
    /* G = Sum (2k+1)!/x^(2k+2), t_0 = 1/x^2. */
    mpfr_ui_div(tg, 1, x2, SRND);
    mpfr_set(g, tg, SRND);
    mpfr_abs(prev, tg, SRND);
    for (unsigned long k = 1; k <= 100000; k++) {
        mpfr_mul_si(tg, tg, (long)((2 * k + 1) * (2 * k)), SRND);
        mpfr_div(tg, tg, x2, SRND);
        mpfr_abs(mag, tg, SRND);
        if (mpfr_cmp(mag, prev) >= 0) break;
        mpfr_add(g, g, tg, SRND);
        mpfr_set(prev, mag, SRND);
        if (mpfr_cmp(mag, thr) < 0) break;
    }

    mpfr_sinh_cosh(sx, cx, x, SRND);
    mpfr_mul(cx, cx, f, SRND);                         /* cosh(x) F */
    mpfr_mul(sx, sx, g, SRND);                         /* sinh(x) G */
    mpfr_add(out, cx, sx, SRND);
    mpfr_clears(f, g, tf, tg, x2, mag, prev, thr, cx, sx, (mpfr_ptr)0);
}

/* Shi(x) for a numeric real `arg` at out_prec bits. Uses odd symmetry to reduce
 * to x >= 0. */
static Expr* shi_mpfr_real(const Expr* arg, mpfr_prec_t out_prec) {
    mpfr_t x;
    mpfr_init2(x, out_prec + 16);
    if (!shi_set_mpfr(x, arg)) { mpfr_clear(x); return NULL; }
    int sgn = mpfr_sgn(x);
    if (sgn == 0) {
        mpfr_clear(x);
        return out_prec <= 53 ? expr_new_real(0.0) : expr_new_mpfr_bits(out_prec);
    }
    double xabs = fabs(mpfr_get_d(x, SRND));

    mpfr_t val;
    if (xabs > (double)(out_prec + 24) * M_LN2 + 8.0) {
        /* asymptotic regime */
        mpfr_prec_t wp = out_prec + 64;
        mpfr_t ax;
        mpfr_inits2(wp, ax, val, (mpfr_ptr)0);
        mpfr_abs(ax, x, SRND);
        shi_real_asymptotic(val, ax, wp);
        mpfr_clear(ax);
    } else {
        /* convergent regime: all-positive terms, no cancellation guard needed */
        mpfr_prec_t wp = out_prec + 64;
        mpfr_t ax;
        mpfr_inits2(wp, ax, val, (mpfr_ptr)0);
        mpfr_abs(ax, x, SRND);
        shi_real_series(val, ax, wp, xabs);
        mpfr_clear(ax);
    }
    if (sgn < 0) mpfr_neg(val, val, SRND);

    Expr* out = shi_emit_component(val, out_prec);
    mpfr_clear(val);
    mpfr_clear(x);
    return out;
}

/* -------------------- complex paths (ncpx) -------------------- */

/* Multiply the components of `z` in place by the integer m. */
static void shi_ncpx_mul_si(ncpx* z, long m) {
    mpfr_mul_si(z->re, z->re, m, SRND);
    mpfr_mul_si(z->im, z->im, m, SRND);
}

/* cosh(w) and sinh(w) into cw, sw via the exponential (ncpx has no hyperbolic
 * primitives): cosh = (e^w + e^-w)/2, sinh = (e^w - e^-w)/2. */
static void shi_ncpx_cosh_sinh(ncpx* cw, ncpx* sw, const ncpx* w, mpfr_prec_t wp) {
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

/* Convergent series Shi(w) = Sum w^(2k+1)/((2k+1)(2k+1)!) into `out`, for w with
 * Re(w) >= 0. All terms positive. Returns false if it fails to converge. */
static bool shi_ncpx_series(ncpx* out, const ncpx* w, mpfr_prec_t wp, double wabs) {
    ncpx term, w2, sum;
    ncpx_init(&term, wp); ncpx_init(&w2, wp); ncpx_init(&sum, wp);
    ncpx_set(&term, w);                    /* k = 0 term = w */
    ncpx_set(&sum, w);
    ncpx_mul(&w2, w, w, wp);

    mpfr_t mag, peak, thr;
    mpfr_inits2(wp, mag, peak, thr, (mpfr_ptr)0);
    ncpx_abs(peak, w);
    long shift = (long)(wp > 8 ? wp - 8 : 1);
    unsigned long cap = (unsigned long)(2.0 * wabs) + (unsigned long)wp + 1000;
    bool ok = false;
    for (unsigned long k = 1; k <= cap; k++) {
        ncpx_mul(&term, &term, &w2, wp);
        shi_ncpx_mul_si(&term, (long)(2 * k - 1));
        mpfr_div_si(term.re, term.re, (long)(2 * k), SRND);
        mpfr_div_si(term.im, term.im, (long)(2 * k), SRND);
        mpfr_div_si(term.re, term.re, (long)(2 * k + 1), SRND);
        mpfr_div_si(term.im, term.im, (long)(2 * k + 1), SRND);
        mpfr_div_si(term.re, term.re, (long)(2 * k + 1), SRND);
        mpfr_div_si(term.im, term.im, (long)(2 * k + 1), SRND);
        ncpx_add(&sum, &sum, &term);
        ncpx_abs(mag, &term);
        if (mpfr_cmp(mag, peak) > 0) mpfr_set(peak, mag, SRND);
        if ((double)(2 * k) > wabs) {
            mpfr_mul_2si(thr, peak, -shift, SRND);
            if (mpfr_cmp(mag, thr) < 0) { ok = true; break; }
        }
    }
    if (ok) ncpx_set(out, &sum);
    mpfr_clears(mag, peak, thr, (mpfr_ptr)0);
    ncpx_clear(&term); ncpx_clear(&w2); ncpx_clear(&sum);
    return ok;
}

/* Sum one auxiliary asymptotic series into `acc` (already = t_0), advancing the
 * running term `t` by the factor (a)(b)/w2 each step to its smallest term.
 * `odd` selects the F-series (base 2k, 2k-1) vs G-series (2k+1, 2k). All terms
 * positive (no sign alternation). */
static void shi_ncpx_aux(ncpx* acc, ncpx* t, const ncpx* w2, mpfr_prec_t wp,
                         mpfr_t thr, int odd) {
    mpfr_t mag, prev;
    mpfr_inits2(wp, mag, prev, (mpfr_ptr)0);
    ncpx_abs(prev, t);
    for (unsigned long k = 1; k <= 100000; k++) {
        long a = odd ? (long)(2 * k) : (long)(2 * k + 1);
        long b = odd ? (long)(2 * k - 1) : (long)(2 * k);
        shi_ncpx_mul_si(t, a);
        shi_ncpx_mul_si(t, b);
        ncpx_div(t, t, w2, wp);
        ncpx_abs(mag, t);
        if (mpfr_cmp(mag, prev) >= 0) break;
        ncpx_add(acc, acc, t);
        mpfr_set(prev, mag, SRND);
        if (mpfr_cmp(mag, thr) < 0) break;
    }
    mpfr_clears(mag, prev, (mpfr_ptr)0);
}

/* Bare asymptotic B(w) = cosh(w) F(w) + sinh(w) G(w) into `out`, for Re(w) >= 0.
 * The caller applies the odd fold and adds the Stokes constant. */
static void shi_ncpx_asymptotic(ncpx* out, const ncpx* w, mpfr_prec_t wp) {
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
    shi_ncpx_aux(&f, &tf, &w2, wp, thr, 1);
    ncpx_div(&tg, &one, &w2, wp);      /* 1/w^2  */
    ncpx_set(&g, &tg);
    shi_ncpx_aux(&g, &tg, &w2, wp, thr, 0);

    shi_ncpx_cosh_sinh(&cw, &sw, w, wp);
    ncpx_mul(&cw, &cw, &f, wp);        /* cosh(w) F */
    ncpx_mul(&sw, &sw, &g, wp);        /* sinh(w) G */
    ncpx_add(out, &cw, &sw);

    mpfr_clear(thr);
    ncpx_clear(&w2); ncpx_clear(&one); ncpx_clear(&f); ncpx_clear(&g);
    ncpx_clear(&tf); ncpx_clear(&tg); ncpx_clear(&cw); ncpx_clear(&sw);
}

/* Build a complex result from (re, im) at out_prec, overflow-aware per component.
 * make_complex drops a zero imaginary part to a bare real. */
static Expr* shi_complex_result(const mpfr_t re, const mpfr_t im, mpfr_prec_t out_prec) {
    return make_complex(shi_emit_component(re, out_prec),
                        shi_emit_component(im, out_prec));
}

/* Shi(z) for a numeric complex `arg` at out_prec bits. Uses odd symmetry to
 * reduce to Re >= 0, then picks the convergent or asymptotic path and adds the
 * i (Pi/2) sign(Im z) Stokes constant on the asymptotic branch. */
static Expr* shi_mpfr_complex(Expr* re, Expr* im, mpfr_prec_t out_prec) {
    double red = 0.0, imd = 0.0;
    (void)shi_to_double(re, &red);
    (void)shi_to_double(im, &imd);
    int sgn = (red < 0.0) ? -1 : 1;      /* Shi(-z) = -Shi(z): fold to Re >= 0 */
    double wre = fabs(red), wim = imd;    /* w = sgn * z has Re >= 0 */
    double wabs = sqrt(wre * wre + wim * wim);

    Expr* out = NULL;
    if (wabs > (double)(out_prec + 24) * M_LN2 + 8.0) {
        mpfr_prec_t wp = out_prec + 96;
        ncpx w, g;
        ncpx_init(&w, wp); ncpx_init(&g, wp);
        if (shi_set_mpfr(w.re, re) && shi_set_mpfr(w.im, im)) {
            if (sgn < 0) ncpx_neg(&w, &w);        /* w = sgn*z, Re(w) >= 0 */
            shi_ncpx_asymptotic(&g, &w, wp);      /* B(w) */
            if (sgn < 0) ncpx_neg(&g, &g);        /* B is odd: B(z) = sgn B(w) */
            if (imd != 0.0) {
                mpfr_t hp;
                mpfr_init2(hp, wp);
                mpfr_const_pi(hp, SRND);
                mpfr_div_2ui(hp, hp, 1, SRND);    /* Pi/2 */
                if (imd < 0.0) mpfr_neg(hp, hp, SRND);
                mpfr_add(g.im, g.im, hp, SRND);   /* + i (Pi/2) sign(Im z) */
                mpfr_clear(hp);
            }
            if (!mpfr_nan_p(g.re) && !mpfr_nan_p(g.im) &&
                !mpfr_inf_p(g.re) && !mpfr_inf_p(g.im))
                out = shi_complex_result(g.re, g.im, out_prec);
        }
        ncpx_clear(&w); ncpx_clear(&g);
        return out;
    }

    /* Convergent series: partial sums peak ~e^|z|, answer ~e^|Im z|; add
     * (|z| + |Re z|)/ln2 guard bits to absorb the cancellation exactly. */
    long guard = 64 + (long)((wabs + wre) / M_LN2);
    mpfr_prec_t wp = out_prec + (mpfr_prec_t)guard;
    ncpx w, g;
    ncpx_init(&w, wp); ncpx_init(&g, wp);
    if (shi_set_mpfr(w.re, re) && shi_set_mpfr(w.im, im)) {
        if (sgn < 0) ncpx_neg(&w, &w);
        if (shi_ncpx_series(&g, &w, wp, wabs)) {
            if (sgn < 0) ncpx_neg(&g, &g);
            if (!mpfr_nan_p(g.re) && !mpfr_nan_p(g.im) &&
                !mpfr_inf_p(g.re) && !mpfr_inf_p(g.im))
                out = shi_complex_result(g.re, g.im, out_prec);
        }
    }
    ncpx_clear(&w); ncpx_clear(&g);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* SinhIntegral[z]                                                    */
/* ------------------------------------------------------------------ */

static Expr* shi_one_arg(Expr* arg) {
    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return expr_new_integer(0);                          /* Shi[0] = 0 */
    if (shi_is_symbol(arg, "Infinity"))         return shi_make_signed_infinity(1);
    if (shi_is_neg_infinity(arg))               return shi_make_signed_infinity(-1);
    if (shi_is_symbol(arg, "ComplexInfinity"))  return expr_new_symbol(SYM_Indeterminate);
    if (shi_is_symbol(arg, "Indeterminate"))    return expr_new_symbol(SYM_Indeterminate);
    {
        int s = shi_directed_imag_infinity(arg);
        if (s != 0) return shi_make_directed_imag_half_pi(s);  /* +-I Pi/2 */
    }

    /* 2. Numeric real (machine and arbitrary precision). */
    if (arg->type == EXPR_REAL) {
#ifdef USE_MPFR
        return shi_mpfr_real(arg, 53);
#else
        double x = arg->data.real, v;
        if (shi_machine_real(x, &v)) return expr_new_real(v);
        return NULL;
#endif
    }
#ifdef USE_MPFR
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        return shi_mpfr_real(arg, prec);
    }
#endif

    /* 3. Numeric complex (Complex[..] with an inexact part). */
    {
        Expr *re, *im;
        if (is_complex(arg, &re, &im) && (shi_is_inexact(re) || shi_is_inexact(im))) {
#ifdef USE_MPFR
            long bits = numeric_min_inexact_bits(arg);
            mpfr_prec_t out_prec = (bits && bits > 53) ? (mpfr_prec_t)bits : 53;
            Expr* out = shi_mpfr_complex(re, im, out_prec);
            if (out) return out;
#else
            double rr, ii;
            if (shi_to_double(re, &rr) && shi_to_double(im, &ii)) {
                double complex v;
                if (shi_machine_complex(rr + ii * I, &v)) {
                    if (cimag(v) == 0.0) return expr_new_real(creal(v));
                    return make_complex(expr_new_real(creal(v)),
                                        expr_new_real(cimag(v)));
                }
            }
#endif
        }
    }

    /* 4. Symbolic odd symmetry: SinhIntegral[-x] -> -SinhIntegral[x]. */
    if (shi_is_neg_leading_times(arg)) {
        Expr* pos = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1), expr_copy(arg) }, 2));
        Expr* inner = expr_new_function(expr_new_symbol(SYM_SinhIntegral), &pos, 1);
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1), inner }, 2));
    }

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

static Expr* shi_emit_argx(size_t argc) {
    fprintf(stderr,
            "SinhIntegral::argx: SinhIntegral called with %zu arguments; "
            "1 argument is expected.\n",
            argc);
    return NULL;
}

Expr* builtin_sinhintegral(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return shi_one_arg(args[0]);
    return shi_emit_argx(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void sinhintegral_init(void) {
    symtab_add_builtin("SinhIntegral", builtin_sinhintegral);
    symtab_get_def("SinhIntegral")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

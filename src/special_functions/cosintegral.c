/* Mathilda -- the cosine integral  Ci(z) = -Int_z^Inf Cos[t]/t dt.
 *
 *   CosIntegral[z]   Ci(z)
 *
 * Ci has a logarithmic singularity at z = 0 and a branch cut running along
 * the negative real axis (-Infinity, 0]. It is NOT entire and NOT odd, so --
 * unlike SinIntegral -- there is no odd-symmetry fold; instead negative /
 * left-half-plane inputs are handled through the principal Log branch and an
 * explicit reflection. Evaluation is layered so each argument takes the
 * cheapest route:
 *
 *   exact special values     ->  -Infinity, 0, I Pi, Infinity, Indeterminate
 *   machine real             ->  MPFR series/asymptotic at 53 bits
 *   arbitrary real           ->  the same, at the input precision
 *   complex (any precision)  ->  the ncpx series/asymptotic with guard bits
 *   everything else          ->  stays symbolic (return NULL)
 *
 * The convergent Maclaurin series (DLMF 6.6.5) is
 *
 *   Ci(z) = EulerGamma + Log(z) + Sum_{k>=1} (-1)^k z^(2k) / (2k (2k)!),
 *
 * valid on the principal branch for all z != 0 -- the principal Log(z) already
 * supplies the correct +-i Pi jump across the cut, so no folding is needed for
 * the convergent path. Its partial sums can reach magnitude ~e^|z| before the
 * O(1)-sized answer emerges, so the MPFR paths add ~|z|/ln2 guard bits to
 * absorb that cancellation exactly. For large |z| the convergent series is
 * infeasible; there we use the asymptotic expansion (DLMF 6.12.4)
 *
 *   Ci(z) = sin(z) f(z) - cos(z) g(z)   (+ Stokes constant, see below),
 *     f(z) ~ Sum (-1)^k (2k)!   / z^(2k+1) = 1/z - 2!/z^3 + ...,
 *     g(z) ~ Sum (-1)^k (2k+1)! / z^(2k+2) = 1/z^2 - 3!/z^4 + ...,
 *
 * (the same f, g SinIntegral computes; only the combination differs -- Si uses
 * Pi/2 - cos f - sin g). The bare sin f - cos g asymptotic is valid in the open
 * right half plane; on/left of the imaginary axis a piecewise-constant Stokes
 * term restores the principal branch:
 *
 *   Ci(z) = sin(z) f(z) - cos(z) g(z) + i C,
 *     C =  0                for Re(z) > 0,
 *     C = +-Pi/2            for Re(z) = 0 (sign of Im z),
 *     C = +-Pi              for Re(z) < 0 (sign of Im z; +Pi from above on the cut).
 *
 * (Equivalently, for Re(z) < 0 we fold w = -z into the right half plane, use
 * the reflection Ci(z) = Ci(-z) + i Pi sign(Im z), and evaluate Ci(-z) by the
 * clean right-half asymptotic.) For a negative real x the result is the
 * from-above limit Complex[Ci(|x|), Pi], matching Ci(-Infinity) = I Pi.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "cosintegral.h"
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
#define CI_EULER_GAMMA 0.57721566490153286061

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

static bool ci_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

static bool ci_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

static bool ci_to_double(const Expr* e, double* out) {
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
static bool ci_is_neg_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!ci_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == -1 && ci_is_symbol(b, "Infinity"))
        return true;
    if (b->type == EXPR_INTEGER && b->data.integer == -1 && ci_is_symbol(a, "Infinity"))
        return true;
    return false;
}

/* Match Times[Complex[0, s], Infinity] (s = +-1): +-I Infinity. Returns
 * +1 for I Infinity, -1 for -I Infinity, 0 if not of this form. */
static int ci_directed_imag_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return 0;
    if (!ci_is_symbol(e->data.function.head, "Times")) return 0;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    Expr *cplx = NULL, *inf = NULL;
    if (ci_is_symbol(b, "Infinity")) { cplx = a; inf = b; }
    else if (ci_is_symbol(a, "Infinity")) { cplx = b; inf = a; }
    if (!inf) return 0;
    Expr *re, *im;
    if (!is_complex(cplx, &re, &im)) return 0;
    if (!(re->type == EXPR_INTEGER && re->data.integer == 0)) return 0;
    if (im->type == EXPR_INTEGER && im->data.integer == 1)  return 1;
    if (im->type == EXPR_INTEGER && im->data.integer == -1) return -1;
    return 0;
}

/* Build -Infinity = Times[-1, Infinity]. */
static Expr* ci_make_neg_infinity(void) {
    return expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2);
}

/* Build I Pi = Times[Complex[0, 1], Pi]. */
static Expr* ci_make_i_pi(void) {
    Expr* imag = make_complex(expr_new_integer(0), expr_new_integer(1));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ imag, expr_new_symbol(SYM_Pi) }, 2));
}

#ifndef USE_MPFR
/* ------------------------------------------------------------------ */
/* Machine-precision fallback paths (USE_MPFR=0 builds)               */
/* ------------------------------------------------------------------ */

/* Ci(a) for a machine real a > 0 via convergent series / asymptotic. */
static bool ci_machine_pos(double a, double* out) {
    double v;
    if (a <= 40.0) {
        double a2 = a * a, term = -a2 / 4.0, sum = -a2 / 4.0;   /* k = 1 term */
        for (int k = 2; k <= 100000; k++) {
            term *= -a2 * (double)(k - 1)
                    / ((double)k * (2 * k) * (2 * k - 1));
            sum += term;
            if ((double)(2 * k) > a && fabs(term) <= 1e-18 * (fabs(sum) + 1.0)) break;
        }
        v = CI_EULER_GAMMA + log(a) + sum;
    } else {
        double f = 1.0 / a, g = 1.0 / (a * a), a2 = a * a;
        double tf = f, tg = g, pf = f, pg = g;
        for (int k = 1; k <= (int)a; k++) {
            tf *= -(double)(2 * k) * (2 * k - 1) / a2;
            if (fabs(tf) >= fabs(pf)) break;
            f += tf; pf = tf;
        }
        for (int k = 1; k <= (int)a; k++) {
            tg *= -(double)(2 * k + 1) * (2 * k) / a2;
            if (fabs(tg) >= fabs(pg)) break;
            g += tg; pg = tg;
        }
        v = sin(a) * f - cos(a) * g;
    }
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

/* Ci(x) for a machine real x != 0. For x < 0 the from-above branch value is
 * Ci(|x|) + i Pi, so (*im) carries Pi. */
static bool ci_machine_real(double x, double* re, double* im) {
    double v;
    if (!ci_machine_pos(fabs(x), &v)) return false;
    *re = v;
    *im = (x < 0.0) ? M_PI : 0.0;
    return true;
}

/* Ci(z) for a machine double-complex z via the convergent series. */
static bool ci_machine_complex(double complex z, double complex* out) {
    double complex z2 = z * z, term = -z2 / 4.0, sum = -z2 / 4.0;   /* k = 1 term */
    double zabs = cabs(z);
    for (int k = 2; k <= 200000; k++) {
        term *= -z2 * (double)(k - 1) / ((double)k * (2 * k) * (2 * k - 1));
        sum += term;
        if ((double)(2 * k) > zabs && cabs(term) <= 1e-17 * (cabs(sum) + 1.0)) break;
    }
    double complex v = CI_EULER_GAMMA + clog(z) + sum;   /* clog: principal */
    if (!isfinite(creal(v)) || !isfinite(cimag(v))) return false;
    *out = v;
    return true;
}
#endif /* !USE_MPFR */

#ifdef USE_MPFR
/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool ci_set_mpfr(mpfr_t out, const Expr* e) {
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

/* -------------------- real paths -------------------- */

/* Convergent even series Sum_{k>=1} (-1)^k a^(2k)/(2k (2k)!) into `sum`
 * (precision wp), for a > 0. Stops once past the term peak (2k > a) and the
 * term has decayed below 2^-(wp-8) of the largest term seen. */
static void ci_real_series(mpfr_t sum, const mpfr_t a, mpfr_prec_t wp, double aabs) {
    mpfr_t term, a2, mag, peak, thr;
    mpfr_inits2(wp, term, a2, mag, peak, thr, (mpfr_ptr)0);
    mpfr_mul(a2, a, a, SRND);
    mpfr_div_ui(term, a2, 4, SRND);       /* k = 1 term = -a^2/4 */
    mpfr_neg(term, term, SRND);
    mpfr_set(sum, term, SRND);
    mpfr_abs(peak, term, SRND);
    long shift = (long)(wp > 8 ? wp - 8 : 1);
    unsigned long cap = (unsigned long)(2.0 * aabs) + (unsigned long)wp + 1000;
    for (unsigned long k = 2; k <= cap; k++) {
        /* term *= -a^2 (k-1) / (k (2k)(2k-1)). */
        mpfr_mul(term, term, a2, SRND);
        mpfr_neg(term, term, SRND);
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

/* Asymptotic Ci(a) = sin(a) f(a) - cos(a) g(a) into `out` (precision wp),
 * for a > 0. Each of f, g is summed to its smallest term. */
static void ci_real_asymptotic(mpfr_t out, const mpfr_t a, mpfr_prec_t wp) {
    mpfr_t f, g, tf, tg, a2, mag, prev, thr, ca, sa;
    mpfr_inits2(wp, f, g, tf, tg, a2, mag, prev, thr, ca, sa, (mpfr_ptr)0);
    mpfr_mul(a2, a, a, SRND);
    mpfr_set_ui(thr, 1, SRND);
    mpfr_mul_2si(thr, thr, -(long)wp, SRND);          /* 2^-wp */

    /* f = Sum (-1)^k (2k)!/a^(2k+1), t_0 = 1/a. */
    mpfr_ui_div(tf, 1, a, SRND);
    mpfr_set(f, tf, SRND);
    mpfr_abs(prev, tf, SRND);
    for (unsigned long k = 1; k <= 100000; k++) {
        mpfr_mul_si(tf, tf, (long)((2 * k) * (2 * k - 1)), SRND);
        mpfr_div(tf, tf, a2, SRND);
        mpfr_neg(tf, tf, SRND);
        mpfr_abs(mag, tf, SRND);
        if (mpfr_cmp(mag, prev) >= 0) break;          /* terms growing */
        mpfr_add(f, f, tf, SRND);
        mpfr_set(prev, mag, SRND);
        if (mpfr_cmp(mag, thr) < 0) break;
    }
    /* g = Sum (-1)^k (2k+1)!/a^(2k+2), t_0 = 1/a^2. */
    mpfr_ui_div(tg, 1, a2, SRND);
    mpfr_set(g, tg, SRND);
    mpfr_abs(prev, tg, SRND);
    for (unsigned long k = 1; k <= 100000; k++) {
        mpfr_mul_si(tg, tg, (long)((2 * k + 1) * (2 * k)), SRND);
        mpfr_div(tg, tg, a2, SRND);
        mpfr_neg(tg, tg, SRND);
        mpfr_abs(mag, tg, SRND);
        if (mpfr_cmp(mag, prev) >= 0) break;
        mpfr_add(g, g, tg, SRND);
        mpfr_set(prev, mag, SRND);
        if (mpfr_cmp(mag, thr) < 0) break;
    }

    mpfr_sin_cos(sa, ca, a, SRND);
    mpfr_mul(sa, sa, f, SRND);                         /* sin(a) f */
    mpfr_mul(ca, ca, g, SRND);                         /* cos(a) g */
    mpfr_sub(out, sa, ca, SRND);                       /* sin f - cos g */
    mpfr_clears(f, g, tf, tg, a2, mag, prev, thr, ca, sa, (mpfr_ptr)0);
}

/* Ci(a) for a > 0 into already-init'd `val` (precision >= out_prec + guard). */
static void ci_real_core(mpfr_t val, const mpfr_t a, double aabs, mpfr_prec_t out_prec) {
    if (aabs > (double)(out_prec + 24) * M_LN2 + 8.0) {
        ci_real_asymptotic(val, a, mpfr_get_prec(val));
    } else {
        mpfr_t euler, lna, s;
        mpfr_prec_t wp = mpfr_get_prec(val);
        mpfr_inits2(wp, euler, lna, s, (mpfr_ptr)0);
        ci_real_series(s, a, wp, aabs);
        mpfr_const_euler(euler, SRND);
        mpfr_log(lna, a, SRND);
        mpfr_add(val, euler, lna, SRND);
        mpfr_add(val, val, s, SRND);                   /* gamma + ln a + series */
        mpfr_clears(euler, lna, s, (mpfr_ptr)0);
    }
}

/* Build a real-or-complex result (re, im) at out_prec: <= 53 yields Real parts,
 * higher yields MPFR. make_complex drops a zero imaginary part to a bare real. */
static Expr* ci_complex_result(const mpfr_t re, const mpfr_t im, mpfr_prec_t out_prec) {
    Expr *rr, *ii;
    if (out_prec <= 53) {
        rr = expr_new_real(mpfr_get_d(re, SRND));
        ii = expr_new_real(mpfr_get_d(im, SRND));
    } else {
        rr = expr_new_mpfr_bits(out_prec);
        ii = expr_new_mpfr_bits(out_prec);
        mpfr_set(rr->data.mpfr, re, SRND);
        mpfr_set(ii->data.mpfr, im, SRND);
    }
    return make_complex(rr, ii);
}

/* Ci(x) for a numeric real `arg` at out_prec bits. x > 0 yields a bare real;
 * x < 0 yields the from-above branch value Complex[Ci(|x|), Pi]; x == 0 yields
 * -Infinity. */
static Expr* ci_mpfr_real(const Expr* arg, mpfr_prec_t out_prec) {
    mpfr_t x;
    mpfr_init2(x, out_prec + 16);
    if (!ci_set_mpfr(x, arg)) { mpfr_clear(x); return NULL; }
    int sgn = mpfr_sgn(x);
    if (sgn == 0) { mpfr_clear(x); return ci_make_neg_infinity(); }

    double aabs = fabs(mpfr_get_d(x, SRND));
    long guard = 64 + (long)(aabs / M_LN2);
    if (aabs > (double)(out_prec + 24) * M_LN2 + 8.0) guard = 64;
    mpfr_prec_t wp = out_prec + (mpfr_prec_t)guard;

    mpfr_t a, val;
    mpfr_inits2(wp, a, val, (mpfr_ptr)0);
    mpfr_abs(a, x, SRND);
    ci_real_core(val, a, aabs, out_prec);

    Expr* out;
    if (sgn > 0) {
        if (out_prec <= 53) out = expr_new_real(mpfr_get_d(val, SRND));
        else { out = expr_new_mpfr_bits(out_prec); mpfr_set(out->data.mpfr, val, SRND); }
    } else {
        mpfr_t pi;
        mpfr_init2(pi, wp);
        mpfr_const_pi(pi, SRND);                        /* from-above jump: +Pi */
        out = ci_complex_result(val, pi, out_prec);
        mpfr_clear(pi);
    }
    mpfr_clears(a, val, (mpfr_ptr)0);
    mpfr_clear(x);
    return out;
}

/* -------------------- complex paths (ncpx) -------------------- */

/* Multiply the components of `z` in place by the integer m. */
static void ci_ncpx_mul_si(ncpx* z, long m) {
    mpfr_mul_si(z->re, z->re, m, SRND);
    mpfr_mul_si(z->im, z->im, m, SRND);
}

/* Divide the components of `z` in place by the positive integer m. */
static void ci_ncpx_div_ui(ncpx* z, unsigned long m) {
    mpfr_div_ui(z->re, z->re, m, SRND);
    mpfr_div_ui(z->im, z->im, m, SRND);
}

/* Convergent series Ci(z) = gamma + Log(z) + Sum (-1)^k z^(2k)/(2k (2k)!) into
 * `out`, valid on the principal branch for any z != 0. Returns false if the
 * even series fails to converge within the cap. */
static bool ci_ncpx_series(ncpx* out, const ncpx* z, mpfr_prec_t wp, double zabs) {
    ncpx term, z2, sum, lg;
    ncpx_init(&term, wp); ncpx_init(&z2, wp);
    ncpx_init(&sum, wp);  ncpx_init(&lg, wp);
    ncpx_mul(&z2, z, z, wp);
    ncpx_set(&term, &z2);                  /* k = 1 term = -z^2/4 */
    ci_ncpx_div_ui(&term, 4);
    ncpx_neg(&term, &term);
    ncpx_set(&sum, &term);

    mpfr_t mag, peak, thr, euler;
    mpfr_inits2(wp, mag, peak, thr, euler, (mpfr_ptr)0);
    ncpx_abs(peak, &term);
    long shift = (long)(wp > 8 ? wp - 8 : 1);
    unsigned long cap = (unsigned long)(2.0 * zabs) + (unsigned long)wp + 1000;
    bool ok = false;
    for (unsigned long k = 2; k <= cap; k++) {
        /* term *= -z^2 (k-1) / (k (2k)(2k-1)). */
        ncpx_mul(&term, &term, &z2, wp);
        ncpx_neg(&term, &term);
        ci_ncpx_mul_si(&term, (long)(k - 1));
        ci_ncpx_div_ui(&term, k);
        ci_ncpx_div_ui(&term, 2 * k);
        ci_ncpx_div_ui(&term, 2 * k - 1);
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
 * running term `t` by the factor -(a)(b)/w2 each step to its smallest term.
 * `odd` selects the f-series (base 2k, 2k-1) vs g-series (2k+1, 2k). */
static void ci_ncpx_aux(ncpx* acc, ncpx* t, const ncpx* w2, mpfr_prec_t wp,
                        mpfr_t thr, int odd) {
    mpfr_t mag, prev;
    mpfr_inits2(wp, mag, prev, (mpfr_ptr)0);
    ncpx_abs(prev, t);
    for (unsigned long k = 1; k <= 100000; k++) {
        long a = odd ? (long)(2 * k) : (long)(2 * k + 1);
        long b = odd ? (long)(2 * k - 1) : (long)(2 * k);
        ci_ncpx_mul_si(t, a);
        ci_ncpx_mul_si(t, b);
        ncpx_div(t, t, w2, wp);
        ncpx_neg(t, t);
        ncpx_abs(mag, t);
        if (mpfr_cmp(mag, prev) >= 0) break;
        ncpx_add(acc, acc, t);
        mpfr_set(prev, mag, SRND);
        if (mpfr_cmp(mag, thr) < 0) break;
    }
    mpfr_clears(mag, prev, (mpfr_ptr)0);
}

/* Bare asymptotic Ci(w) = sin(w) f(w) - cos(w) g(w) into `out`, for Re(w) >= 0.
 * The caller adds the Stokes constant. */
static void ci_ncpx_asymptotic(ncpx* out, const ncpx* w, mpfr_prec_t wp) {
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
    ci_ncpx_aux(&f, &tf, &w2, wp, thr, 1);
    ncpx_div(&tg, &one, &w2, wp);      /* 1/w^2  */
    ncpx_set(&g, &tg);
    ci_ncpx_aux(&g, &tg, &w2, wp, thr, 0);

    ncpx_sin(&sw, w, wp);
    ncpx_cos(&cw, w, wp);
    ncpx_mul(&sw, &sw, &f, wp);        /* sin(w) f */
    ncpx_mul(&cw, &cw, &g, wp);        /* cos(w) g */
    ncpx_sub(out, &sw, &cw);           /* sin f - cos g */

    mpfr_clear(thr);
    ncpx_clear(&w2); ncpx_clear(&one); ncpx_clear(&f); ncpx_clear(&g);
    ncpx_clear(&tf); ncpx_clear(&tg); ncpx_clear(&cw); ncpx_clear(&sw);
}

/* Ci(z) for a numeric complex `arg` at out_prec bits. Moderate |z| uses the
 * principal convergent series (correct across the whole plane); large |z| uses
 * the right-half asymptotic with the Stokes / reflection constant. */
static Expr* ci_mpfr_complex(Expr* re, Expr* im, mpfr_prec_t out_prec) {
    double red = 0.0, imd = 0.0;
    (void)ci_to_double(re, &red);
    (void)ci_to_double(im, &imd);
    double zabs = sqrt(red * red + imd * imd);

    Expr* out = NULL;
    if (zabs > (double)(out_prec + 24) * M_LN2 + 8.0) {
        /* Asymptotic. Fold Re(z) < 0 to the right half via the reflection
         * Ci(z) = Ci(-z) + i Pi sign(Im z); add the Stokes constant otherwise. */
        mpfr_prec_t wp = out_prec + 96;
        int fold = (red < 0.0);
        ncpx w, g;
        ncpx_init(&w, wp); ncpx_init(&g, wp);
        if (ci_set_mpfr(w.re, re) && ci_set_mpfr(w.im, im)) {
            if (fold) ncpx_neg(&w, &w);                 /* w = -z, Re(w) > 0 */
            ci_ncpx_asymptotic(&g, &w, wp);
            /* Stokes / reflection constant added to the imaginary part. */
            if (red != 0.0 || imd != 0.0) {
                mpfr_t pi;
                mpfr_init2(pi, wp);
                mpfr_const_pi(pi, SRND);
                if (red > 0.0) {
                    /* C = 0 */
                } else if (red == 0.0) {
                    mpfr_div_2ui(pi, pi, 1, SRND);       /* Pi/2 */
                    if (imd < 0.0) mpfr_neg(pi, pi, SRND);
                    mpfr_add(g.im, g.im, pi, SRND);
                } else {
                    if (imd < 0.0) mpfr_neg(pi, pi, SRND);
                    mpfr_add(g.im, g.im, pi, SRND);      /* +- Pi */
                }
                mpfr_clear(pi);
            }
            if (!mpfr_nan_p(g.re) && !mpfr_nan_p(g.im) &&
                !mpfr_inf_p(g.re) && !mpfr_inf_p(g.im))
                out = ci_complex_result(g.re, g.im, out_prec);
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
    if (ci_set_mpfr(z.re, re) && ci_set_mpfr(z.im, im)) {
        if (ci_ncpx_series(&g, &z, wp, zabs)) {
            if (!mpfr_nan_p(g.re) && !mpfr_nan_p(g.im) &&
                !mpfr_inf_p(g.re) && !mpfr_inf_p(g.im))
                out = ci_complex_result(g.re, g.im, out_prec);
        }
    }
    ncpx_clear(&z); ncpx_clear(&g);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* CosIntegral[z]                                                     */
/* ------------------------------------------------------------------ */

static Expr* ci_one_arg(Expr* arg) {
    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return ci_make_neg_infinity();                       /* Ci[0] = -Inf */
    if (ci_is_symbol(arg, "Infinity"))         return expr_new_integer(0);   /* Ci[Inf]=0 */
    if (ci_is_neg_infinity(arg))               return ci_make_i_pi();        /* Ci[-Inf]=I Pi */
    if (ci_is_symbol(arg, "ComplexInfinity"))  return expr_new_symbol(SYM_Indeterminate);
    if (ci_is_symbol(arg, "Indeterminate"))    return expr_new_symbol(SYM_Indeterminate);
    {
        int s = ci_directed_imag_infinity(arg);
        if (s != 0) return expr_new_symbol(SYM_Infinity);    /* Ci[+-I Inf] = Infinity */
    }

    /* 2. Numeric real (machine and arbitrary precision). */
    if (arg->type == EXPR_REAL) {
#ifdef USE_MPFR
        return ci_mpfr_real(arg, 53);
#else
        double x = arg->data.real;
        if (x == 0.0) return ci_make_neg_infinity();
        double rr, ii;
        if (ci_machine_real(x, &rr, &ii)) {
            if (ii == 0.0) return expr_new_real(rr);
            return make_complex(expr_new_real(rr), expr_new_real(ii));
        }
        return NULL;
#endif
    }
#ifdef USE_MPFR
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        return ci_mpfr_real(arg, prec);
    }
#endif

    /* 3. Numeric complex (Complex[..] with an inexact part). */
    {
        Expr *re, *im;
        if (is_complex(arg, &re, &im) && (ci_is_inexact(re) || ci_is_inexact(im))) {
#ifdef USE_MPFR
            long bits = numeric_min_inexact_bits(arg);
            mpfr_prec_t out_prec = (bits && bits > 53) ? (mpfr_prec_t)bits : 53;
            Expr* out = ci_mpfr_complex(re, im, out_prec);
            if (out) return out;
#else
            double rr, ii;
            if (ci_to_double(re, &rr) && ci_to_double(im, &ii)) {
                double complex v;
                if (ci_machine_complex(rr + ii * I, &v)) {
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

static Expr* ci_emit_argx(size_t argc) {
    fprintf(stderr,
            "CosIntegral::argx: CosIntegral called with %zu arguments; "
            "1 argument is expected.\n",
            argc);
    return NULL;
}

Expr* builtin_cosintegral(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return ci_one_arg(args[0]);
    return ci_emit_argx(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void cosintegral_init(void) {
    symtab_add_builtin("CosIntegral", builtin_cosintegral);
    symtab_get_def("CosIntegral")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

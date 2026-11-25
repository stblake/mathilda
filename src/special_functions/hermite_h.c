/* Mathilda -- HermiteH, the physicists' Hermite polynomial H_n(x).
 *
 *   HermiteH[n, x]
 *
 * The Hermite polynomials satisfy  y'' - 2x y' + 2n y = 0  and are orthogonal on
 * (-Inf, Inf) with weight e^(-x^2). Evaluation is layered so each argument shape
 * takes the cheapest exact or numeric route (mirroring the ChebyshevU design in
 * chebyshev_u.c, whose integer-coefficient, two-argument substrate HermiteH
 * shares via chebyshev_common.{c,h}):
 *
 *   non-negative integer order n
 *     symbolic / exact-number x  ->  the explicit degree-n polynomial in x,
 *                                    with exact integer coefficients, from
 *                                    FLINT's fmpz_poly_hermite_h (or, without
 *                                    FLINT, the three-term recurrence); an exact
 *                                    numeric x then evaluates to an exact number.
 *     inexact x                  ->  the stable three-term recurrence evaluated
 *                                    in machine/MPFR complex arithmetic, NOT by
 *                                    substituting into the polynomial (whose
 *                                    cancellation ruins accuracy for large n,
 *                                    e.g. HermiteH[100, N[5,16]]).
 *   non-integer / complex / negative order, some inexact argument
 *                                ->  H_n(z) = 2^n Sqrt[Pi] (
 *                                        1F1(-n/2; 1/2; z^2) / Gamma[(1-n)/2]
 *                                      - 2z 1F1((1-n)/2; 3/2; z^2) / Gamma[-n/2]),
 *                                    built as an expression over
 *                                    Hypergeometric1F1 so it inherits its
 *                                    machine + MPFR, real + complex numerics.
 *   special value
 *     H_0(x) = 1,  H_n(0) = 2^n Sqrt[Pi] / Gamma[(1-n)/2].
 *   everything else              ->  stays symbolic (return NULL).
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * The machine-precision double / double-complex recurrence lives in sf_machine.c
 * (sf_machine_hermite_h[_c]) and is shared bit-for-bit with the NDArray
 * element-wise kernel and the Compile[] VM (ndkernels.c). The scalar machine
 * path here dispatches to the same helpers.
 *
 * Series / SeriesCoefficient: the derivative rule D[HermiteH[n,z],z] =
 * 2 n HermiteH[n-1,z] lives in calculus/deriv.c and, with H_n(0) above, gives
 * the Taylor expansion at a finite point for free via series.c's Taylor-via-D;
 * the asymptotic expansion at Infinity and the symbolic-index SeriesCoefficient
 * are bespoke kernels in calculus/series.c.
 *
 * Deferred (left symbolic): negative-integer order with an exact argument (no
 * Erfc closed form; a numeric argument is handled by the 1F1 form above), and
 * bespoke tight Interval enclosures for non-integer order (the generic
 * NumericFunction + D[] interval threading in eval.c handles those best-effort;
 * an integer order threads exactly through the explicit polynomial).
 */
#include "hermite_h.h"
#include "chebyshev_common.h"   /* cheb_* shared substrate (family-agnostic) */

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include <gmp.h>

#include "arithmetic.h"      /* make_rational, make_complex */
#include "attr.h"
#include "common.h"          /* builtin_arg_error */
#include "eval.h"            /* eval_and_free */
#include "expr.h"
#include "numeric.h"         /* get_approx_mpfr */
#include "numeric_complex.h" /* ncpx toolkit */
#include "sf_machine.h"      /* sf_machine_hermite_h[_c] -- raw-double path */
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#define PRND MPFR_RNDN
#endif

#ifdef USE_FLINT
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#endif

/* Largest order for which the explicit polynomial is built; beyond this a
 * symbolic/exact-x call stays symbolic. */
#define HERMITE_POLY_CAP 2000

/* ------------------------------------------------------------------ */
/* Expression-building shorthands                                      */
/* ------------------------------------------------------------------ */

static Expr* mk_fn1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* mk_fn3(const char* h, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b, c }, 3);
}

/* ------------------------------------------------------------------ */
/* Explicit integer-coefficient polynomial                             */
/* ------------------------------------------------------------------ */

/* Fill `out` (length deg+1, mpz_init'd to 0) with the integer coefficients of
 * H_deg(x) = Sum_i out[i] x^i. */
static void hermite_coeffs(unsigned long deg, mpz_t* out) {
#ifdef USE_FLINT
    fmpz_poly_t p;
    fmpz_poly_init(p);
    fmpz_poly_hermite_h(p, deg);
    fmpz_t fc;
    fmpz_init(fc);
    for (unsigned long i = 0; i <= deg; i++) {
        fmpz_poly_get_coeff_fmpz(fc, p, (slong)i);
        fmpz_get_mpz(out[i], fc);
    }
    fmpz_clear(fc);
    fmpz_poly_clear(p);
#else
    /* Three-term recurrence H_0 = 1, H_1 = 2x, H_k = 2x H_{k-1} - 2(k-1) H_{k-2}:
     * new[i] = 2 cur[i-1] - 2(k-1) prev[i]. */
    if (deg == 0) { mpz_set_ui(out[0], 1); return; }
    mpz_t* prev = cheb_mpz_array_new(deg + 1);   /* H_{k-2} */
    mpz_t* cur  = cheb_mpz_array_new(deg + 1);   /* H_{k-1} */
    mpz_t* nw   = cheb_mpz_array_new(deg + 1);   /* H_k     */
    if (!prev || !cur || !nw) {
        cheb_mpz_array_free(prev, deg + 1);
        cheb_mpz_array_free(cur, deg + 1);
        cheb_mpz_array_free(nw, deg + 1);
        return;
    }
    mpz_set_ui(prev[0], 1);                       /* H_0 = 1  */
    mpz_set_ui(cur[1], 2);                        /* H_1 = 2x */
    if (deg == 1) {
        for (size_t i = 0; i <= deg; i++) mpz_set(out[i], cur[i]);
    } else {
        for (unsigned long k = 2; k <= deg; k++) {
            unsigned long c2 = 2UL * (k - 1);
            for (unsigned long i = 0; i <= deg; i++) {
                mpz_set_ui(nw[i], 0);
                if (i >= 1) mpz_mul_2exp(nw[i], cur[i - 1], 1);  /* 2 cur[i-1]     */
                mpz_submul_ui(nw[i], prev[i], c2);               /* - 2(k-1) prev[i] */
            }
            for (unsigned long i = 0; i <= deg; i++) {
                mpz_set(prev[i], cur[i]);
                mpz_set(cur[i], nw[i]);
            }
        }
        for (size_t i = 0; i <= deg; i++) mpz_set(out[i], cur[i]);
    }
    cheb_mpz_array_free(prev, deg + 1);
    cheb_mpz_array_free(cur, deg + 1);
    cheb_mpz_array_free(nw, deg + 1);
#endif
}

/* The explicit polynomial H_deg(x). */
static Expr* hermite_poly(unsigned long deg, const Expr* x) {
    mpz_t* c = cheb_mpz_array_new(deg + 1);
    if (!c) return NULL;
    hermite_coeffs(deg, c);
    Expr* poly = cheb_poly_from_coeffs(c, deg, x);
    cheb_mpz_array_free(c, deg + 1);
    return poly;
}

/* ------------------------------------------------------------------ */
/* Special value and closed forms                                      */
/* ------------------------------------------------------------------ */

/* H_n(0) = 2^n Sqrt[Pi] / Gamma[(1-n)/2]. Built as an expression and evaluated,
 * so it serves a symbolic n (-> the Gamma form) and an inexact n (-> a number)
 * alike. */
static Expr* hermite_value_at_zero(const Expr* n_e) {
    Expr* pow2 = mk_fn2("Power", expr_new_integer(2), expr_copy((Expr*)n_e));
    Expr* sqrtpi = mk_fn2("Power", expr_new_symbol("Pi"), make_rational(1, 2));
    Expr* one_minus_n = mk_fn2("Plus", expr_new_integer(1),
                           mk_fn2("Times", expr_new_integer(-1), expr_copy((Expr*)n_e)));
    Expr* half_arg = mk_fn2("Times", make_rational(1, 2), one_minus_n);  /* (1-n)/2 */
    Expr* invgam = mk_fn2("Power", mk_fn1("Gamma", half_arg), expr_new_integer(-1));
    Expr* res = expr_new_function(expr_new_symbol("Times"),
                   (Expr*[]){ pow2, sqrtpi, invgam }, 3);
    return eval_and_free(res);
}

/* Build the closed form
 *   H_n(z) = 2^n Sqrt[Pi] ( 1F1(-n/2; 1/2; z^2)/Gamma[(1-n)/2]
 *                         - 2z 1F1((1-n)/2; 3/2; z^2)/Gamma[-n/2] )
 * as an unevaluated expression over Hypergeometric1F1 (HermiteH is entire in z,
 * so one expression covers the whole plane). The caller evaluates it; see
 * hermite_numeric_hyp for the precision handling. */
static Expr* hermite_hyp_form(const Expr* n_e, const Expr* x_e) {
    Expr* z2 = mk_fn2("Power", expr_copy((Expr*)x_e), expr_new_integer(2));   /* z^2 */
    /* term1 = 1F1(-n/2; 1/2; z^2) / Gamma[(1-n)/2] */
    Expr* a1 = mk_fn2("Times", make_rational(-1, 2), expr_copy((Expr*)n_e));  /* -n/2 */
    Expr* hyp1 = mk_fn3("Hypergeometric1F1", a1, make_rational(1, 2), expr_copy(z2));
    Expr* half1mn = mk_fn2("Times", make_rational(1, 2),
                        mk_fn2("Plus", expr_new_integer(1),
                            mk_fn2("Times", expr_new_integer(-1), expr_copy((Expr*)n_e))));
    Expr* gam1 = mk_fn1("Gamma", expr_copy(half1mn));                         /* Gamma[(1-n)/2] */
    Expr* term1 = mk_fn2("Times", hyp1, mk_fn2("Power", gam1, expr_new_integer(-1)));
    /* term2 = -2z 1F1((1-n)/2; 3/2; z^2) / Gamma[-n/2] */
    Expr* hyp2 = mk_fn3("Hypergeometric1F1", half1mn, make_rational(3, 2), z2);
    Expr* nn2 = mk_fn2("Times", make_rational(-1, 2), expr_copy((Expr*)n_e)); /* -n/2 */
    Expr* gam2 = mk_fn1("Gamma", nn2);                                        /* Gamma[-n/2] */
    Expr* term2 = expr_new_function(expr_new_symbol("Times"),
                     (Expr*[]){ expr_new_integer(-2), expr_copy((Expr*)x_e), hyp2,
                                mk_fn2("Power", gam2, expr_new_integer(-1)) }, 4);
    Expr* inner = mk_fn2("Plus", term1, term2);
    return expr_new_function(expr_new_symbol("Times"),
               (Expr*[]){ mk_fn2("Power", expr_new_integer(2), expr_copy((Expr*)n_e)),
                          mk_fn2("Power", expr_new_symbol("Pi"), make_rational(1, 2)),
                          inner }, 3);
}

#ifdef USE_MPFR
/* Extra working precision (bits) to absorb the catastrophic cancellation in the
 * 1F1-difference form: the two confluent-hypergeometric terms each grow like
 * e^(Re z^2) while H_n(z) itself is O(|2z|^Re n), so ~Re(z^2) log2(e) bits
 * cancel. A fixed guard (as elsewhere) is not enough for a large real z (e.g.
 * HermiteH[1/2, 5.] loses ~16 digits), so the guard scales with Re(z^2) (times
 * two for the series' own intermediate overshoot), floored and capped. */
static long hermite_hyp_guard_bits(const Expr* x_e) {
    double xr = 0.0, xi = 0.0;
    mpfr_t re, im;
    mpfr_init2(re, 64); mpfr_init2(im, 64);
    bool inexact = false;
    if (get_approx_mpfr(x_e, re, im, &inexact)) {
        xr = mpfr_get_d(re, PRND);
        xi = mpfr_get_d(im, PRND);
    }
    mpfr_clear(re); mpfr_clear(im);
    double re_z2 = xr * xr - xi * xi;             /* Re(z^2) */
    if (!(re_z2 > 0.0)) re_z2 = 0.0;              /* also catches NaN */
    long g = (long)(re_z2 * 1.4426950408889634 * 2.0) + 64;
    if (g > 8192) g = 8192;                       /* bound the cost for huge |z| */
    return g;
}
#endif

/* Numeric H_n(z) for a non-integer / complex / negative order via the closed
 * form, built as an expression so it inherits Hypergeometric1F1's machine /
 * MPFR, real / complex numerics. With MPFR, the form is evaluated at a working
 * precision raised by a cancellation guard (hermite_hyp_guard_bits) and then
 * rounded back to the precision the arguments imply, so a large real argument
 * does not silently lose digits. Returns NULL if it did not collapse to a
 * number. */
/* True when x is a real number > 0 (machine, MPFR or exact). */
static bool hermite_real_positive(const Expr* x_e) {
#ifdef USE_MPFR
    mpfr_t re, im;
    mpfr_init2(re, 64); mpfr_init2(im, 64);
    bool inexact = false, ok = get_approx_mpfr(x_e, re, im, &inexact)
                              && mpfr_zero_p(im) && mpfr_sgn(re) > 0;
    mpfr_clear(re); mpfr_clear(im);
    return ok;
#else
    /* no MPFR: a machine real or an exact positive integer / rational */
    int64_t n, d;
    if (x_e->type == EXPR_REAL) return x_e->data.real > 0.0;
    if (x_e->type == EXPR_INTEGER) return x_e->data.integer > 0;
    if (is_rational((Expr*)x_e, &n, &d)) return n > 0;
    return false;
#endif
}

/* H_n(x) = 2^n U(-n/2, 1/2, x^2) for Re x > 0 (DLMF 12.7.14 with 13.2.42-type
 * relation; Tricomi U). Unlike the 1F1 difference it has no cancellation: the
 * difference loses ~2 Sqrt(|n| x^2) log2(e) bits for a large negative order,
 * which the Re(x^2)-scaled guard does not cover -- HermiteH[-100., 10.] had
 * six correct digits, HermiteH[-1000., 10.] came out -1.0*^-1358 where the
 * value is +2.57*^-1608, and HermiteH[-1.*^10, 10.] an exact 0 (numeric-tower
 * review G10, checked against the integral representation). */
static Expr* hermite_u_form(const Expr* n_e, const Expr* x_e) {
    Expr* a = mk_fn2("Times", make_rational(-1, 2), expr_copy((Expr*)n_e));   /* -n/2 */
    Expr* z2 = mk_fn2("Power", expr_copy((Expr*)x_e), expr_new_integer(2));   /* x^2 */
    Expr* u = mk_fn3("HypergeometricU", a, make_rational(1, 2), z2);
    return mk_fn2("Times", mk_fn2("Power", expr_new_integer(2), expr_copy((Expr*)n_e)), u);
}

static Expr* hermite_numeric_hyp_form(const Expr* n_e, const Expr* x_e, bool use_u);

#ifndef USE_MPFR
/* H_{-m}(x) for a real order -m < 0 and a real x, in a build without MPFR:
 * there HypergeometricU has no numerics, so the U form above cannot apply and
 * the 1F1 difference cancels to nothing (HermiteH[-100., 10.] was 0. where the
 * value is 1.31*^-138). The Laplace-type integral representation, valid for
 * Re nu < 0 and every x,
 *     H_{-m}(x) = 1/Gamma(m) Integral_0^Inf e^(-t^2 - 2 x t) t^(m-1) dt
 * (checked against the MPFR build's U form, x of either sign, m from 1e-3 to
 * 1e6), has a positive integrand, so nothing cancels. With t = t0 e^y about its peak
 * t0 = (-x + Sqrt[x^2 + 2m])/2 the integrand is entire in y and decays on both
 * sides, so the trapezoid rule on the whole line converges geometrically: the
 * step min(sigma/2, 1/16) (sigma the peak width) keeps the discretisation error
 * below 1e-17. Far left, where t^2 + 2|x| t is below an ulp, the terms are
 * exactly geometric and their tail is summed in closed form (expm1, not
 * 1 - exp, or a tiny m loses everything). The value is assembled in logs; its
 * accuracy is that of the problem's conditioning in m (~1e-14 at m = 100).
 * False (decline) outside the domain or on a non-finite result. */
static bool hermite_neg_order_machine(double m, double x, double* out) {
    if (!(m > 0.0) || !isfinite(m) || !isfinite(x)) return false;
    double t0 = (x >= 0.0) ? m / (x + sqrt(x * x + 2.0 * m))     /* no cancellation */
                           : 0.5 * (-x + sqrt(x * x + 2.0 * m));
    double sigma = 1.0 / sqrt(2.0 * t0 * t0 + m);
    double h = fmin(sigma / 2.0, 0.0625);
    double s = 1.0, c = 0.0;                                    /* y = 0 term */
    long steps = 0;
    for (int dir = -1; dir <= 1; dir += 2) {
        for (long k = 1; ; k++) {
            if (++steps > 4000000L) return false;
            double y = (double)dir * (double)k * h, ey = exp(y), t = t0 * ey;
            double g = -t0 * t0 * (ey * ey - 1.0) - 2.0 * x * t0 * (ey - 1.0) + m * y;
            double term = exp(g);
            double s2 = s + term, bb = s2 - s;                  /* compensated */
            c += (s - (s2 - bb)) + (term - bb); s = s2;
            if (g < -50.0 && term < 1e-20 * s) break;
            if (dir < 0 && t * t + 2.0 * fabs(x) * t < 1e-18) {  /* geometric from here */
                double tail = term * exp(-m * h) / -expm1(-m * h);
                s2 = s + tail; bb = s2 - s; c += (s - (s2 - bb)) + (tail - bb); s = s2;
                break;
            }
        }
    }
    s += c;
    double v = exp(m * log(t0) - t0 * t0 - 2.0 * x * t0 + log(s * h) - lgamma(m));
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

/* A real machine order < 0 (machine Real, or an exact negative number). */
static bool hermite_neg_order_value(const Expr* n_e, double* m) {
    int64_t p, q;
    if (n_e->type == EXPR_REAL) *m = -n_e->data.real;
    else if (n_e->type == EXPR_INTEGER) *m = -(double)n_e->data.integer;
    else if (is_rational((Expr*)n_e, &p, &q)) *m = -(double)p / (double)q;
    else return false;
    return *m > 0.0;
}
#endif

static Expr* hermite_numeric_hyp(const Expr* n_e, const Expr* x_e) {
#ifndef USE_MPFR
    {
        double m, xr, xi, v; bool real;
        if (hermite_neg_order_value(n_e, &m) && cheb_machine_xy(x_e, &xr, &xi, &real) && real
            && hermite_neg_order_machine(m, xr, &v))
            return expr_new_real(v);
    }
#endif
    /* The U form where it applies; the 1F1 form when it does not collapse to
     * a number (a build whose HypergeometricU numerics are absent: without
     * FLINT, HermiteH[-100., 10.] stayed unevaluated). */
    if (hermite_real_positive(x_e)) {
        Expr* r = hermite_numeric_hyp_form(n_e, x_e, true);
        if (r) return r;
    }
    return hermite_numeric_hyp_form(n_e, x_e, false);
}

static Expr* hermite_numeric_hyp_form(const Expr* n_e, const Expr* x_e, bool use_u) {
    Expr* cf = use_u ? hermite_u_form(n_e, x_e) : hermite_hyp_form(n_e, x_e);
#ifdef USE_MPFR
    long op = cheb_out_prec(n_e, x_e);            /* target output precision, >=53 */
    NumericSpec hi;
    hi.mode = NUMERIC_MODE_MPFR;
    /* the U form needs no cancellation guard, only the usual few bits */
    hi.bits = op + (use_u ? 32 : hermite_hyp_guard_bits(x_e));
    Expr* hiout = numericalize(cf, hi);           /* evaluate with guard digits */
    expr_free(cf);
    if (!hiout || !cheb_is_number_leaf(hiout)) { expr_free(hiout); return NULL; }
    /* Round to the precision the inputs imply: a machine double when they are
     * machine, else an op-bit MPFR value (capping, never inventing, digits). */
    NumericSpec tgt;
    if (op <= 53) {
        tgt = numeric_machine_spec();
    } else {
        tgt.mode = NUMERIC_MODE_MPFR_CAP;
        tgt.bits = op;
    }
    Expr* out = numericalize(hiout, tgt);
    expr_free(hiout);
    if (out && cheb_is_number_leaf(out)) return out;
    expr_free(out);
    return NULL;
#else
    Expr* out = eval_and_free(cf);
    if (out && cheb_is_number_leaf(out)) return out;
    expr_free(out);
    return NULL;
#endif
}

/* ------------------------------------------------------------------ */
/* Numeric evaluation (MPFR recurrence)                                */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* H_deg(x) for an inexact x, by the stable forward recurrence in complex MPFR
 * arithmetic. Returns NULL if x is not numeric. */
static Expr* hermite_numeric_recurrence(unsigned long deg, const Expr* x_e,
                                        mpfr_prec_t out_prec) {
    mpfr_prec_t wp = (out_prec < 64 ? 64 : out_prec) + 32;

    ncpx x, h0, h1, hk, t, s, c;
    ncpx_init(&x, wp);  ncpx_init(&h0, wp); ncpx_init(&h1, wp);
    ncpx_init(&hk, wp); ncpx_init(&t, wp);  ncpx_init(&s, wp); ncpx_init(&c, wp);

    Expr* out = NULL;
    bool inexact = false;
    if (get_approx_mpfr(x_e, x.re, x.im, &inexact)) {
        bool real_input = mpfr_zero_p(x.im);

        ncpx_set_ui(&h0, 1);        /* H_0 = 1  */
        ncpx_add(&h1, &x, &x);      /* H_1 = 2x */

        ncpx* res;
        if      (deg == 0) res = &h0;
        else if (deg == 1) res = &h1;
        else {
            for (unsigned long k = 2; k <= deg; k++) {
                ncpx_mul(&t, &x, &h1, wp);        /* t = x h1            */
                ncpx_add(&t, &t, &t);             /* t = 2x h1           */
                ncpx_set_ui(&c, 2UL * (k - 1));   /* c = 2(k-1)          */
                ncpx_mul(&s, &c, &h0, wp);        /* s = 2(k-1) h0       */
                ncpx_sub(&hk, &t, &s);            /* hk = 2x h1-2(k-1)h0 */
                ncpx_set(&h0, &h1);
                ncpx_set(&h1, &hk);
            }
            res = &h1;
        }
        if (!mpfr_nan_p(res->re) && !mpfr_nan_p(res->im))
            out = cheb_make_result(res->re, res->im, out_prec, real_input);
    }

    ncpx_clear(&x);  ncpx_clear(&h0); ncpx_clear(&h1);
    ncpx_clear(&hk); ncpx_clear(&t);  ncpx_clear(&s); ncpx_clear(&c);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Two-argument driver                                                 */
/* ------------------------------------------------------------------ */

static Expr* hermite_h_two_arg(Expr* n_e, Expr* x_e) {
    int64_t ni;
    bool n_is_int = cheb_exact_int(n_e, &ni);
    bool n_nonneg_int = n_is_int && ni >= 0;
    bool x_inexact = cheb_is_inexact(x_e);

    /* Special value at x = 0 for a non-(non-negative-integer) order. A
     * non-negative integer order is handled exactly by the polynomial below
     * (HermiteH[1,0] = 0, HermiteH[0,0] = 1, ...). Holds for symbolic and
     * inexact n, since the built expression evaluates either way. */
    /* H_n(0) = 0 for every odd n, a bignum included (the Gamma form below
     * overflowed for HermiteH[2^70 + 1, 0], HermiteH#2). */
    if (n_e->type == EXPR_BIGINT && mpz_sgn(n_e->data.bigint) > 0 &&
        mpz_odd_p(n_e->data.bigint) && cheb_is_int(x_e, 0))
        return expr_new_integer(0);
    if (!n_nonneg_int && cheb_is_int(x_e, 0))
        return hermite_value_at_zero(n_e);

#ifdef USE_MPFR
    /* Numeric non-negative integer order with an inexact argument: the stable
     * recurrence. Covers an exact integer order and an inexact integer-valued
     * one alike. */
    if (x_inexact) {
        int64_t k = -1;
        bool have = n_nonneg_int ? (k = ni, true)
                                 : (cheb_inexact_int_order(n_e, &k) && k >= 0);
        if (have && k >= 0 && k <= (int64_t)CHEB_NUMERIC_DEG_CAP) {
            unsigned long deg = (unsigned long)k;
            mpfr_prec_t op = cheb_out_prec(n_e, x_e);
            /* Machine fast path: the raw-double / double-complex recurrence,
             * bit-identical to the ND kernel and Compile[]. Taken only at
             * machine precision (<= 53 bits) and degree, with a machine
             * real/complex x; a large degree, overflow, or arbitrary-precision
             * argument makes the helper decline, and we fall through to the MPFR
             * recurrence below unchanged. */
            if (op <= 53 && deg <= 1000) {
                double xr, xi; bool real;
                if (cheb_machine_xy(x_e, &xr, &xi, &real)) {
                    if (real) {
                        double v;
                        if (sf_machine_hermite_h((double)deg, xr, &v))
                            return expr_new_real(v);
                    } else {
                        double rr, ri;
                        if (sf_machine_hermite_h_c((double)deg, xr, xi, &rr, &ri))
                            return make_complex(expr_new_real(rr), expr_new_real(ri));
                    }
                }
            }
            Expr* out = hermite_numeric_recurrence(deg, x_e, op);
            if (out) return out;
        }
    }
#endif

    /* Exact non-negative integer order with a symbolic / exact-number x: the
     * explicit polynomial (exact for an exact x). Too large stays symbolic. */
    if (n_nonneg_int) {
        if ((unsigned long)ni > HERMITE_POLY_CAP) return NULL;
        return hermite_poly((unsigned long)ni, x_e);
    }

    /* Non-integer / complex / negative order: numeric only when an argument is
     * inexact. */
    if (cheb_is_inexact(n_e) || x_inexact) {
        Expr* out = hermite_numeric_hyp(n_e, x_e);
        if (out) return out;
    }
    return NULL;                                     /* stay symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                 */
/* ------------------------------------------------------------------ */

Expr* builtin_hermite_h(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return builtin_arg_error("HermiteH", argc, 2, 2);
    Expr** args = res->data.function.args;
    return hermite_h_two_arg(args[0], args[1]);
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

void hermite_h_init(void) {
    symtab_add_builtin("HermiteH", builtin_hermite_h);
    symtab_get_def("HermiteH")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

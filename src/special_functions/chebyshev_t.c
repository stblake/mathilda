/* Mathilda -- ChebyshevT, the Chebyshev polynomial of the first kind T_n(x).
 *
 *   ChebyshevT[n, x]
 *
 * Evaluation mirrors ChebyshevU (chebyshev_u.c): each argument shape takes the
 * cheapest exact or numeric route.
 *
 *   exact integer order n
 *     symbolic / exact-number x  ->  the explicit degree-|n| polynomial in x,
 *                                    with exact integer coefficients, from
 *                                    FLINT's fmpz_poly_chebyshev_t (or, without
 *                                    FLINT, the three-term recurrence); an exact
 *                                    numeric x then evaluates to an exact number.
 *     inexact x                  ->  the stable three-term recurrence evaluated
 *                                    in machine/MPFR complex arithmetic, NOT by
 *                                    substituting into the polynomial (whose
 *                                    cancellation ruins accuracy for large n --
 *                                    e.g. ChebyshevT[60, N[2/3,16]]).
 *   non-integer / complex order, some inexact argument
 *                                ->  T_v(x) = 2F1(-v, v; 1/2; (1-x)/2), built as
 *                                    an expression over Hypergeometric2F1 so it
 *                                    inherits its machine + MPFR, real + complex
 *                                    numerics; for |(1-x)/2| >= 1 the whole-plane
 *                                    form T_v(x) = Cos[v ArcCos[x]] takes over.
 *   special values
 *     T_0(x) = 1,  T_n(0) = Cos[n Pi/2],  T_n(1) = 1.
 *   everything else              ->  stays symbolic (return NULL).
 *
 * Unlike U, T is even in its order index -- T_{-n}(x) = T_n(x) -- so a negative
 * order reflects to |n| with no sign and no special zero case (T_{-1}(x) = x).
 *
 * The recurrence seed is the only arithmetic difference from U: T_1 = x where
 * U_1 = 2x; the step T_k = 2x T_{k-1} - T_{k-2} is identical.
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Deferred (left symbolic), matching ChebyshevU: symbolic Series /
 * SeriesCoefficient, asymptotic expansions, half-integer closed forms, exact
 * Interval[{-1,1}] at infinite order, and the Cosh[..]-style symbolic
 * simplification of a fully-exact complex-order call (e.g.
 * ChebyshevT[2-7I, 5+3I] stays symbolic). D[ChebyshevT[n,x],x] = n U_{n-1}(x)
 * lives in calculus/deriv.c.
 */
#include "chebyshev_t.h"
#include "chebyshev_common.h"
#include "sym_names.h"

#include <limits.h>
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
#include "sf_machine.h"      /* sf_machine_chebyshev_t[_c] -- raw-double path */
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#define PRND MPFR_RNDN
#endif

#ifdef USE_FLINT
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#endif

/* ------------------------------------------------------------------ */
/* Explicit integer-coefficient polynomial                             */
/* ------------------------------------------------------------------ */

/* Fill `out` (length deg+1, mpz_init'd to 0) with the integer coefficients of
 * T_deg(x) = Sum_i out[i] x^i. */
static void cheb_coeffs(unsigned long deg, mpz_t* out) {
#ifdef USE_FLINT
    fmpz_poly_t p;
    fmpz_poly_init(p);
    fmpz_poly_chebyshev_t(p, deg);
    fmpz_t fc;
    fmpz_init(fc);
    for (unsigned long i = 0; i <= deg; i++) {
        fmpz_poly_get_coeff_fmpz(fc, p, (slong)i);
        fmpz_get_mpz(out[i], fc);
    }
    fmpz_clear(fc);
    fmpz_poly_clear(p);
#else
    /* Three-term recurrence T_0 = 1, T_1 = x, T_k = 2x T_{k-1} - T_{k-2}:
     * new[i] = 2 cur[i-1] - prev[i]  (same step as U; only the T_1 seed differs). */
    if (deg == 0) { mpz_set_ui(out[0], 1); return; }
    mpz_t* prev = cheb_mpz_array_new(deg + 1);   /* T_{k-2} */
    mpz_t* cur  = cheb_mpz_array_new(deg + 1);   /* T_{k-1} */
    mpz_t* nw   = cheb_mpz_array_new(deg + 1);   /* T_k     */
    if (!prev || !cur || !nw) {
        cheb_mpz_array_free(prev, deg + 1);
        cheb_mpz_array_free(cur, deg + 1);
        cheb_mpz_array_free(nw, deg + 1);
        return;
    }
    mpz_set_ui(prev[0], 1);                       /* T_0 = 1 */
    mpz_set_ui(cur[1], 1);                        /* T_1 = x */
    if (deg == 1) {
        for (size_t i = 0; i <= deg; i++) mpz_set(out[i], cur[i]);
    } else {
        for (unsigned long k = 2; k <= deg; k++) {
            for (unsigned long i = 0; i <= deg; i++) {
                mpz_set_ui(nw[i], 0);
                if (i >= 1) mpz_mul_2exp(nw[i], cur[i - 1], 1);  /* 2 cur[i-1] */
                mpz_sub(nw[i], nw[i], prev[i]);                  /* - prev[i]  */
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

/* The explicit polynomial T_deg(x). No sign: T is even in its order index, so
 * the negative-order reflection is handled entirely by |n| in cheb_reflect. */
static Expr* cheb_poly(unsigned long deg, const Expr* x) {
    mpz_t* c = cheb_mpz_array_new(deg + 1);
    if (!c) return NULL;
    cheb_coeffs(deg, c);
    Expr* poly = cheb_poly_from_coeffs(c, deg, x);
    cheb_mpz_array_free(c, deg + 1);
    return poly;
}

/* Reflect an integer order to a non-negative degree: T_{-n}(x) = T_n(x), so the
 * degree is simply |ni| and there is no sign. Returns false for an order whose
 * |value| exceeds the numeric cap (checked before negation, so no overflow). */
static bool cheb_reflect(int64_t ni, unsigned long* deg) {
    if (ni > (int64_t)CHEB_NUMERIC_DEG_CAP || ni < -(int64_t)CHEB_NUMERIC_DEG_CAP)
        return false;
    *deg = (unsigned long)(ni < 0 ? -ni : ni);
    return true;
}

/* ------------------------------------------------------------------ */
/* Numeric evaluation                                                  */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* Above this degree a real arbitrary-precision x takes the closed form
 * cos(n acos x) / cosh(n acosh x) rather than the recurrence. */
#define CHEB_TRIG_MIN_DEG 64

/* T_deg(x) for an inexact x, by the stable forward recurrence (seed T_0 = 1,
 * T_1 = x). A real x runs a real MPFR recurrence; a complex x the complex ncpx
 * one. Splitting the real case is not cosmetic: the complex recurrence carries
 * an always-zero imaginary part through every ncpx_mul (4 real mults/step), so
 * a real arbitrary-precision argument -- the common case, and every deg>1000
 * real argument that falls off the machine fast path -- paid ~2x for nothing.
 * Returns NULL if x is not numeric. */
static Expr* cheb_numeric_recurrence(unsigned long deg,
                                     const Expr* x_e, mpfr_prec_t out_prec) {
    mpfr_prec_t wp = (out_prec < 64 ? 64 : out_prec) + 32;

    ncpx x;
    ncpx_init(&x, wp);

    Expr* out = NULL;
    bool inexact = false;
    if (get_approx_mpfr(x_e, x.re, x.im, &inexact)) {
        if (mpfr_zero_p(x.im) && deg > CHEB_TRIG_MIN_DEG) {
            /* real x, high degree: T_n(x) = cos(n acos x) on [-1, 1], and
             * sgn(x)^n cosh(n acosh |x|) off it -- O(1) operations instead of
             * the recurrence's O(n), and free of its error growth (a rounding
             * at step k reaches T_n multiplied by U_{n-k}(x) ~ n - k near
             * x = +-1, ~n^2/2 ulps in all; the 32 guard bits are gone by
             * n ~ 10^5).  The product n theta gains log2(n) bits of absolute
             * error, which the extra guard bits absorb. */
            mpfr_prec_t gp = wp;
            for (unsigned long m = deg; m > 1; m >>= 1) gp++;
            mpfr_t th, v;
            mpfr_init2(th, gp); mpfr_init2(v, gp);
            mpfr_abs(v, x.re, PRND);
            bool inside = mpfr_cmp_ui(v, 1) <= 0;
            if (inside) {
                mpfr_acos(th, x.re, PRND);
                mpfr_mul_ui(th, th, deg, PRND);
                mpfr_cos(v, th, PRND);
            } else {
                mpfr_acosh(th, v, PRND);
                mpfr_mul_ui(th, th, deg, PRND);
                mpfr_cosh(v, th, PRND);
                if (mpfr_sgn(x.re) < 0 && (deg & 1UL)) mpfr_neg(v, v, PRND);
            }
            if (mpfr_number_p(v))
                out = cheb_make_result(v, v, out_prec, true);
            mpfr_clear(th); mpfr_clear(v);
        } else if (mpfr_zero_p(x.im)) {
            /* real recurrence: plain MPFR, half the per-step work */
            mpfr_t r0, r1, rk, t;
            mpfr_init2(r0, wp); mpfr_init2(r1, wp);
            mpfr_init2(rk, wp); mpfr_init2(t, wp);
            mpfr_set_ui(r0, 1, PRND);       /* T_0 = 1 */
            mpfr_set(r1, x.re, PRND);       /* T_1 = x */

            mpfr_ptr res;
            if      (deg == 0) res = r0;
            else if (deg == 1) res = r1;
            else {
                for (unsigned long k = 2; k <= deg; k++) {
                    mpfr_mul(t, x.re, r1, PRND);  /* t = x t1      */
                    mpfr_mul_2ui(t, t, 1, PRND);  /* t = 2x t1     */
                    mpfr_sub(rk, t, r0, PRND);    /* rk = 2x t1-t0 */
                    mpfr_set(r0, r1, PRND);
                    mpfr_set(r1, rk, PRND);
                }
                res = r1;
            }
            if (!mpfr_nan_p(res))
                out = cheb_make_result(res, res, out_prec, true);

            mpfr_clear(r0); mpfr_clear(r1); mpfr_clear(rk); mpfr_clear(t);
        } else if (deg > CHEB_TRIG_MIN_DEG && deg <= (unsigned long)LONG_MAX) {
            /* complex x, high degree: T_n(z) = (w^n + w^-n)/2 with
             * w = z + Sqrt[z - 1] Sqrt[z + 1] (either root of
             * w^2 - 2 z w + 1 = 0 gives the same sum, as w^-1 is the other
             * one), w^n by binary powering: O(log n) complex products where
             * the recurrence takes n, 1.15 ms against ~0.1 ms for
             * ChebyshevT[2000, z] at 30 digits.  The relative error of w is
             * multiplied by n in w^n, so log2(n) guard bits are added, as for
             * the real closed form above; the two terms only cancel where
             * |w| ~ 1 (z near [-1, 1]), and there to an absolute error of
             * n ulps of |w^n| ~ 1, within the recurrence's own bound. */
            mpfr_prec_t gp = wp + 8;
            for (unsigned long m = deg; m > 1; m >>= 1) gp++;
            ncpx zg, a, b, w, wn, iv;
            ncpx_init(&zg, gp); ncpx_init(&a, gp); ncpx_init(&b, gp);
            ncpx_init(&w, gp); ncpx_init(&wn, gp); ncpx_init(&iv, gp);
            ncpx_set(&zg, &x);
            ncpx one; ncpx_init(&one, gp); ncpx_set_ui(&one, 1);
            ncpx_sub(&a, &zg, &one);  ncpx_sqrt(&a, &a, gp);   /* Sqrt[z - 1] */
            ncpx_add(&b, &zg, &one);  ncpx_sqrt(&b, &b, gp);   /* Sqrt[z + 1] */
            ncpx_mul(&w, &a, &b, gp);
            ncpx_add(&w, &w, &zg);
            /* the root of modulus >= 1, so that w^-n is the small term */
            mpfr_t mw; mpfr_init2(mw, gp);
            ncpx_abs(mw, &w);
            if (mpfr_cmp_ui(mw, 1) < 0) {
                ncpx_div(&w, &one, &w, gp);
            }
            mpfr_clear(mw);
            ncpx_pow_si(&wn, &w, (long)deg, gp);
            ncpx_div(&iv, &one, &wn, gp);
            ncpx_add(&wn, &wn, &iv);
            mpfr_div_2ui(wn.re, wn.re, 1, PRND);
            mpfr_div_2ui(wn.im, wn.im, 1, PRND);
            if (mpfr_number_p(wn.re) && mpfr_number_p(wn.im))
                out = cheb_make_result(wn.re, wn.im, out_prec, false);
            ncpx_clear(&zg); ncpx_clear(&a); ncpx_clear(&b); ncpx_clear(&w);
            ncpx_clear(&wn); ncpx_clear(&iv); ncpx_clear(&one);
        }
        if (!out && !mpfr_zero_p(x.im)) {
            /* complex recurrence in ncpx arithmetic */
            ncpx t0, t1, tk, t;
            ncpx_init(&t0, wp); ncpx_init(&t1, wp);
            ncpx_init(&tk, wp); ncpx_init(&t, wp);
            ncpx_set_ui(&t0, 1);        /* T_0 = 1 */
            ncpx_set(&t1, &x);          /* T_1 = x */

            ncpx* res;
            if      (deg == 0) res = &t0;
            else if (deg == 1) res = &t1;
            else {
                for (unsigned long k = 2; k <= deg; k++) {
                    ncpx_mul(&t, &x, &t1, wp);   /* t = x t1      */
                    ncpx_add(&t, &t, &t);        /* t = 2x t1     */
                    ncpx_sub(&tk, &t, &t0);      /* tk = 2x t1-t0 */
                    ncpx_set(&t0, &t1);
                    ncpx_set(&t1, &tk);
                }
                res = &t1;
            }
            if (!mpfr_nan_p(res->re) && !mpfr_nan_p(res->im))
                out = cheb_make_result(res->re, res->im, out_prec, false);

            ncpx_clear(&t0); ncpx_clear(&t1); ncpx_clear(&tk); ncpx_clear(&t);
        }
    }

    ncpx_clear(&x);
    return out;
}
#endif /* USE_MPFR */

/* T_v(x) = 2F1(-v, v; 1/2; (1-x)/2), built as an expression so it inherits
 * Hypergeometric2F1's machine/MPFR/complex numerics. The Gauss series converges
 * for |(1-x)/2| < 1; returns NULL otherwise (or when it did not collapse to a
 * number), handing off to the ArcCos form below. */
static Expr* cheb_numeric_2f1(const Expr* n_e, const Expr* x_e) {
    /* z = (1 - x)/2 */
    Expr* mx = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)x_e) }, 2);
    Expr* one_minus_x = expr_new_function(expr_new_symbol(SYM_Plus),
                   (Expr*[]){ expr_new_integer(1), mx }, 2);
    Expr* z = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ make_rational(1, 2), one_minus_x }, 2);
    /* a = -v, b = v, c = 1/2 */
    Expr* a = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)n_e) }, 2);
    Expr* b = expr_copy((Expr*)n_e);
    Expr* c = make_rational(1, 2);
    Expr* hyp = expr_new_function(expr_new_symbol(SYM_Hypergeometric2F1),
                   (Expr*[]){ a, b, c, z }, 4);
    Expr* out = eval_and_free(hyp);
    if (out && cheb_is_number_leaf(out)) return out;
    expr_free(out);
    return NULL;
}

/* T_v(x) = Cos[v ArcCos[x]], built as an expression so it inherits ArcCos/Cos's
 * numerics AND their branch conventions -- which match Mathematica's cut from
 * -Infinity to -1 for non-integer v. This covers the whole plane, in particular
 * the |(1-x)/2| >= 1 region the 2F1 series cannot reach. Returns NULL if it did
 * not collapse to a number. */
static Expr* cheb_numeric_arccos(const Expr* n_e, const Expr* x_e) {
    Expr* ac = expr_new_function(expr_new_symbol(SYM_ArcCos),
                   (Expr*[]){ expr_copy((Expr*)x_e) }, 1);
    Expr* arg = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ expr_copy((Expr*)n_e), ac }, 2);
    Expr* res = expr_new_function(expr_new_symbol(SYM_Cos), (Expr*[]){ arg }, 1);
    Expr* out = eval_and_free(res);
    if (out && cheb_is_number_leaf(out)) return out;
    expr_free(out);
    return NULL;
}

/* Non-integer / complex order numeric: the convergent 2F1 form first, then the
 * whole-plane ArcCos closed form for |(1-x)/2| >= 1. */
static Expr* cheb_numeric_general(const Expr* n_e, const Expr* x_e) {
    Expr* out = cheb_numeric_2f1(n_e, x_e);
    if (out) return out;
    return cheb_numeric_arccos(n_e, x_e);
}

/* ------------------------------------------------------------------ */
/* Special values and the two-argument driver                          */
/* ------------------------------------------------------------------ */

static Expr* chebyshev_t_two_arg(Expr* n_e, Expr* x_e) {
    int64_t ni;
    bool n_is_int = cheb_exact_int(n_e, &ni);
    bool x_inexact = cheb_is_inexact(x_e);

    /* Special values for a non-integer order (an integer order is handled
     * exactly by the polynomial / recurrence below). Both hold for symbolic
     * and inexact n, since the built expression evaluates either way. */
    if (!n_is_int) {
        if (cheb_is_int(x_e, 0)) return cheb_value_at_zero(n_e);   /* T_n(0)=Cos[n Pi/2] */
        if (cheb_is_int(x_e, 1)) return expr_new_integer(1);        /* T_n(1)=1 */
    }

#ifdef USE_MPFR
    /* Numeric integer order with an inexact argument: the stable recurrence.
     * Covers an exact integer order and an inexact integer-valued one alike
     * (the latter would otherwise pick up a spurious imaginary part from the
     * general ArcCos form). */
    if (x_inexact) {
        int64_t k = 0;
        bool have = n_is_int ? (k = ni, true) : cheb_inexact_int_order(n_e, &k);
        if (have) {
            unsigned long deg;
            if (cheb_reflect(k, &deg)) {   /* T_{-n} = T_n: degree |k|, no sign */
                /* Machine fast path: the raw-double / double-complex recurrence,
                 * bit-identical to the ND kernel and Compile[].  Taken only when
                 * the requested precision is machine (<= 53 bits) and x is a
                 * machine real/complex.  A real x of any degree has a helper
                 * (the recurrence to degree 1000, the closed form Cos[n ArcCos[x]]
                 * / Cosh[n ArcCosh[x]] above); a complex x of degree > 1000, an
                 * overflow, or an arbitrary-precision argument makes it decline,
                 * and we fall through to the MPFR recurrence below unchanged. */
                if (cheb_out_prec(n_e, x_e) <= 53) {
                    double xr, xi; bool real;
                    if (cheb_machine_xy(x_e, &xr, &xi, &real)) {
                        if (real) {
                            double v;
                            if (sf_machine_chebyshev_t((double)deg, xr, &v))
                                return expr_new_real(v);
                        } else {
                            double rr, ri;
                            if (sf_machine_chebyshev_t_c((double)deg, xr, xi, &rr, &ri))
                                return make_complex(expr_new_real(rr),
                                                    expr_new_real(ri));
                        }
                    }
                }
                Expr* out = cheb_numeric_recurrence(deg, x_e,
                                                    cheb_out_prec(n_e, x_e));
                if (out) return out;
            }
        }
    }
#endif

    /* Exact integer order with a symbolic / exact-number x: the explicit
     * polynomial (exact for an exact x). Too large stays symbolic. */
    if (n_is_int) {
        unsigned long deg;
        if (!cheb_reflect(ni, &deg)) return NULL;                 /* absurd order */
        if (deg > CHEB_POLY_CAP) return NULL;
        return cheb_poly(deg, x_e);
    }

    /* Non-integer / complex order: numeric only when an argument is inexact. */
    if (cheb_is_inexact(n_e) || x_inexact) {
        Expr* out = cheb_numeric_general(n_e, x_e);
        if (out) return out;
    }
    return NULL;                                     /* stay symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                 */
/* ------------------------------------------------------------------ */

Expr* builtin_chebyshev_t(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return builtin_arg_error("ChebyshevT", argc, 2, 2);
    Expr** args = res->data.function.args;
    return chebyshev_t_two_arg(args[0], args[1]);
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

void chebyshev_t_init(void) {
    symtab_add_builtin("ChebyshevT", builtin_chebyshev_t);
    symtab_get_def("ChebyshevT")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

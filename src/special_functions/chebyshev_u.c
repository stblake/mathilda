/* Mathilda -- ChebyshevU, the Chebyshev polynomial of the second kind U_n(x).
 *
 *   ChebyshevU[n, x]
 *
 * Evaluation is layered so each argument shape takes the cheapest exact or
 * numeric route (mirroring the LegendreP design in legendre.c):
 *
 *   exact integer order n
 *     symbolic / exact-number x  ->  the explicit degree-|n'| polynomial in x,
 *                                    with exact integer coefficients, from
 *                                    FLINT's fmpz_poly_chebyshev_u (or, without
 *                                    FLINT, the three-term recurrence); an exact
 *                                    numeric x then evaluates to an exact number.
 *     inexact x                  ->  the stable three-term recurrence evaluated
 *                                    in machine/MPFR complex arithmetic, NOT by
 *                                    substituting into the polynomial (whose
 *                                    cancellation ruins accuracy for large n).
 *   non-integer / complex order, some inexact argument
 *                                ->  U_v(x) = (v+1) 2F1(-v, v+2; 3/2; (1-x)/2),
 *                                    built as an expression over Hypergeometric2F1
 *                                    so it inherits its machine + MPFR, real +
 *                                    complex numerics (and analytic continuation).
 *   special values
 *     U_0(x) = 1,  U_n(0) = Cos[n Pi/2],  U_n(1) = n + 1.
 *   everything else              ->  stays symbolic (return NULL).
 *
 * Negative integer order reflects to a non-negative degree:
 *   U_{-1}(x) = 0,  and  U_{-m}(x) = -U_{m-2}(x)  for m >= 2.
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Deferred (left symbolic): symbolic Series / SeriesCoefficient, D[] of a
 * symbolic order (the D[ChebyshevU[n,x],x] recurrence lives in calculus/deriv.c),
 * half-integer closed forms, and bespoke tight Interval enclosures (the generic
 * NumericFunction + D[] interval threading in eval.c handles those best-effort).
 *
 * The argument-shape predicates, the explicit-polynomial assembler, the MPFR
 * result builder, the output-precision rule, and U_n(0) = Cos[n Pi/2] are shared
 * with ChebyshevT via chebyshev_common.{c,h}; only the recurrence seed, the
 * negative-order reflection, and the non-integer-order closed forms live here.
 */
#include "chebyshev_u.h"
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
#include "sf_machine.h"      /* sf_machine_chebyshev_u[_c] -- raw-double path */
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
 * U_deg(x) = Sum_i out[i] x^i. */
static void cheb_coeffs(unsigned long deg, mpz_t* out) {
#ifdef USE_FLINT
    fmpz_poly_t p;
    fmpz_poly_init(p);
    fmpz_poly_chebyshev_u(p, deg);
    fmpz_t fc;
    fmpz_init(fc);
    for (unsigned long i = 0; i <= deg; i++) {
        fmpz_poly_get_coeff_fmpz(fc, p, (slong)i);
        fmpz_get_mpz(out[i], fc);
    }
    fmpz_clear(fc);
    fmpz_poly_clear(p);
#else
    /* Three-term recurrence U_0 = 1, U_1 = 2x, U_k = 2x U_{k-1} - U_{k-2}:
     * new[i] = 2 cur[i-1] - prev[i]. */
    if (deg == 0) { mpz_set_ui(out[0], 1); return; }
    mpz_t* prev = cheb_mpz_array_new(deg + 1);   /* U_{k-2} */
    mpz_t* cur  = cheb_mpz_array_new(deg + 1);   /* U_{k-1} */
    mpz_t* nw   = cheb_mpz_array_new(deg + 1);   /* U_k     */
    if (!prev || !cur || !nw) {
        cheb_mpz_array_free(prev, deg + 1);
        cheb_mpz_array_free(cur, deg + 1);
        cheb_mpz_array_free(nw, deg + 1);
        return;
    }
    mpz_set_ui(prev[0], 1);                       /* U_0 = 1  */
    mpz_set_ui(cur[1], 2);                        /* U_1 = 2x */
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

/* The explicit polynomial sign * U_deg(x). A sign of -1 folds the negative-order
 * reflection into the coefficients so the result prints distributed (1 - 4 x^2
 * rather than -(-1 + 4 x^2)). */
static Expr* cheb_poly(unsigned long deg, int sign, const Expr* x) {
    mpz_t* c = cheb_mpz_array_new(deg + 1);
    if (!c) return NULL;
    cheb_coeffs(deg, c);
    if (sign < 0)
        for (size_t i = 0; i <= deg; i++) mpz_neg(c[i], c[i]);
    Expr* poly = cheb_poly_from_coeffs(c, deg, x);
    cheb_mpz_array_free(c, deg + 1);
    return poly;
}

/* Reflect an integer order to a non-negative degree and a sign: U_{-m}(x) =
 * -U_{m-2}(x) for m >= 2. The caller handles ni == -1 (U_{-1} = 0) first.
 * Returns false for an order whose |value| exceeds the numeric cap. */
static bool cheb_reflect(int64_t ni, unsigned long* deg, int* sign) {
    if (ni < -1) {
        if (ni < -(int64_t)CHEB_NUMERIC_DEG_CAP) return false;
        *deg = (unsigned long)(-ni - 2);
        *sign = -1;
    } else {
        if (ni > (int64_t)CHEB_NUMERIC_DEG_CAP) return false;
        *deg = (unsigned long)ni;
        *sign = 1;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/* Numeric evaluation                                                  */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* U_deg(x) for an inexact x, by the stable forward recurrence in complex MPFR
 * arithmetic. `sign` negates the result (negative-order reflection). Returns
 * NULL if x is not numeric. */
static Expr* cheb_numeric_recurrence(unsigned long deg, int sign,
                                     const Expr* x_e, mpfr_prec_t out_prec) {
    mpfr_prec_t wp = (out_prec < 64 ? 64 : out_prec) + 32;

    ncpx x, u0, u1, uk, t;
    ncpx_init(&x, wp);  ncpx_init(&u0, wp); ncpx_init(&u1, wp);
    ncpx_init(&uk, wp); ncpx_init(&t, wp);

    Expr* out = NULL;
    bool inexact = false;
    if (get_approx_mpfr(x_e, x.re, x.im, &inexact)) {
        bool real_input = mpfr_zero_p(x.im);

        ncpx_set_ui(&u0, 1);        /* U_0 = 1  */
        ncpx_add(&u1, &x, &x);      /* U_1 = 2x */

        ncpx* res;
        if      (deg == 0) res = &u0;
        else if (deg == 1) res = &u1;
        else {
            for (unsigned long k = 2; k <= deg; k++) {
                ncpx_mul(&t, &x, &u1, wp);   /* t = x u1      */
                ncpx_add(&t, &t, &t);        /* t = 2x u1     */
                ncpx_sub(&uk, &t, &u0);      /* uk = 2x u1-u0 */
                ncpx_set(&u0, &u1);
                ncpx_set(&u1, &uk);
            }
            res = &u1;
        }
        if (sign < 0) ncpx_neg(res, res);
        if (!mpfr_nan_p(res->re) && !mpfr_nan_p(res->im))
            out = cheb_make_result(res->re, res->im, out_prec, real_input);
    }

    ncpx_clear(&x);  ncpx_clear(&u0); ncpx_clear(&u1);
    ncpx_clear(&uk); ncpx_clear(&t);
    return out;
}
#endif /* USE_MPFR */

/* U_v(x) = (v+1) 2F1(-v, v+2; 3/2; (1-x)/2), built as an expression so it
 * inherits Hypergeometric2F1's machine/MPFR/complex numerics. The Gauss series
 * converges for |(1-x)/2| < 1; returns NULL otherwise (or when it did not
 * collapse to a number), handing off to the ArcCos form below. */
static Expr* cheb_numeric_2f1(const Expr* n_e, const Expr* x_e) {
    /* z = (1 - x)/2 */
    Expr* mx = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)x_e) }, 2);
    Expr* one_minus_x = expr_new_function(expr_new_symbol(SYM_Plus),
                   (Expr*[]){ expr_new_integer(1), mx }, 2);
    Expr* z = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ make_rational(1, 2), one_minus_x }, 2);
    /* a = -v, b = v + 2, c = 3/2 */
    Expr* a = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)n_e) }, 2);
    Expr* b = expr_new_function(expr_new_symbol(SYM_Plus),
                   (Expr*[]){ expr_copy((Expr*)n_e), expr_new_integer(2) }, 2);
    Expr* c = make_rational(3, 2);
    Expr* hyp = expr_new_function(expr_new_symbol(SYM_Hypergeometric2F1),
                   (Expr*[]){ a, b, c, z }, 4);
    /* (v + 1) * 2F1 */
    Expr* np1 = expr_new_function(expr_new_symbol(SYM_Plus),
                   (Expr*[]){ expr_copy((Expr*)n_e), expr_new_integer(1) }, 2);
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ np1, hyp }, 2);
    Expr* out = eval_and_free(prod);
    if (out && cheb_is_number_leaf(out)) return out;
    expr_free(out);
    return NULL;
}

/* U_v(x) = Sin[(v+1) ArcCos[x]] / Sin[ArcCos[x]], built as an expression so it
 * inherits ArcCos/Sin's numerics AND their branch conventions -- which match
 * Mathematica's cut from -Infinity to -1 for non-integer v. This covers the
 * whole plane, in particular the |(1-x)/2| >= 1 region the 2F1 series cannot
 * reach. Returns NULL if it did not collapse to a number (e.g. x = +-1, where
 * both Sin terms vanish). */
static Expr* cheb_numeric_arccos(const Expr* n_e, const Expr* x_e) {
    Expr* np1 = expr_new_function(expr_new_symbol(SYM_Plus),
                   (Expr*[]){ expr_copy((Expr*)n_e), expr_new_integer(1) }, 2);
    Expr* ac_num = expr_new_function(expr_new_symbol(SYM_ArcCos),
                   (Expr*[]){ expr_copy((Expr*)x_e) }, 1);
    Expr* arg = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ np1, ac_num }, 2);
    Expr* num = expr_new_function(expr_new_symbol(SYM_Sin), (Expr*[]){ arg }, 1);

    Expr* ac_den = expr_new_function(expr_new_symbol(SYM_ArcCos),
                   (Expr*[]){ expr_copy((Expr*)x_e) }, 1);
    Expr* den = expr_new_function(expr_new_symbol(SYM_Sin), (Expr*[]){ ac_den }, 1);
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ den, expr_new_integer(-1) }, 2);
    Expr* res = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ num, inv }, 2);
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

static Expr* chebyshev_u_two_arg(Expr* n_e, Expr* x_e) {
    int64_t ni;
    bool n_is_int = cheb_exact_int(n_e, &ni);
    bool x_inexact = cheb_is_inexact(x_e);

    /* Special values for a non-integer order (an integer order is handled
     * exactly by the polynomial / recurrence below). Both hold for symbolic
     * and inexact n, since the built expression evaluates either way. */
    if (!n_is_int) {
        if (cheb_is_int(x_e, 0)) return cheb_value_at_zero(n_e);          /* U_n(0)=Cos[n Pi/2] */
        /* U_n(-1) = lim Sin[(n+1) t]/Sin[t] at t = Pi: a pole for a real
         * non-integer n (ChebyshevU#2 of the 2026-11 audit -- the inexact
         * form returned the rounding residue 1/Sin[Pi.] ~ -8*10^15, the exact
         * one leaked Power::infy). */
        {
            double nv;
            bool x_m1 = cheb_is_int(x_e, -1)
                || (x_e->type == EXPR_REAL && x_e->data.real == -1.0)
#ifdef USE_MPFR
                || (x_e->type == EXPR_MPFR && mpfr_cmp_si(x_e->data.mpfr, -1) == 0)
#endif
                ;
            bool n_real = n_e->type == EXPR_REAL
                || (n_e->type == EXPR_FUNCTION && is_rational(n_e, NULL, NULL))
#ifdef USE_MPFR
                || n_e->type == EXPR_MPFR
#endif
                ;
            if (x_m1 && n_real) {
                nv = n_e->type == EXPR_REAL ? n_e->data.real : 0.5;
#ifdef USE_MPFR
                if (n_e->type == EXPR_MPFR) nv = mpfr_get_d(n_e->data.mpfr, MPFR_RNDN);
#endif
                if (nv != floor(nv)) return expr_new_symbol(SYM_ComplexInfinity);
            }
        }
        if (cheb_is_int(x_e, 1))                                         /* U_n(1)=n+1 */
            return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                       (Expr*[]){ expr_copy(n_e), expr_new_integer(1) }, 2));
    }

#ifdef USE_MPFR
    /* An exact integer order past the recurrence's cap, or a bignum one, with
     * a real inexact x: the O(1) closed form (U_{-n} = -U_{n-2}). */
    if (x_inexact && (n_e->type == EXPR_BIGINT ||
                      (n_is_int && (ni > (int64_t)CHEB_NUMERIC_DEG_CAP ||
                                    ni < -(int64_t)CHEB_NUMERIC_DEG_CAP)))) {
        mpz_t n;
        mpz_init(n);
        expr_to_mpz(n_e, n);
        bool negate = mpz_sgn(n) < 0;
        if (negate) { mpz_neg(n, n); mpz_sub_ui(n, n, 2); }
        Expr* v = cheb_closed_form_real(1, n, x_e, cheb_out_prec(n_e, x_e));
        mpz_clear(n);
        if (v && negate)
            v = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                  (Expr*[]){ expr_new_integer(-1), v }, 2));
        if (v) return v;
    }

    /* Numeric integer order with an inexact argument: the stable recurrence.
     * Covers an exact integer order and an inexact integer-valued one alike
     * (the latter would otherwise pick up a spurious imaginary part from the
     * general ArcCos form). */
    if (x_inexact) {
        int64_t k = 0;
        bool have = n_is_int ? (k = ni, true) : cheb_inexact_int_order(n_e, &k);
        if (have) {
            if (k == -1) return expr_new_integer(0);                     /* U_{-1}(x)=0 */
            unsigned long deg; int sign;
            if (cheb_reflect(k, &deg, &sign)) {
                /* Machine fast path: the raw-double / double-complex recurrence,
                 * bit-identical to the ND kernel and Compile[].  Taken only when
                 * the requested precision is machine (<= 53 bits) and x is a
                 * machine real/complex.  A large degree, an overflow, or an
                 * arbitrary-precision argument makes the helper decline, and we
                 * fall through to the MPFR recurrence below unchanged. */
                if (cheb_out_prec(n_e, x_e) <= 53 && deg <= 1000) {
                    double xr, xi; bool real;
                    if (cheb_machine_xy(x_e, &xr, &xi, &real)) {
                        if (real) {
                            double v;
                            if (sf_machine_chebyshev_u((double)deg, xr, &v))
                                return expr_new_real(sign < 0 ? -v : v);
                        } else {
                            double rr, ri;
                            if (sf_machine_chebyshev_u_c((double)deg, xr, xi, &rr, &ri))
                                return make_complex(
                                    expr_new_real(sign < 0 ? -rr : rr),
                                    expr_new_real(sign < 0 ? -ri : ri));
                        }
                    }
                }
                Expr* out = cheb_numeric_recurrence(deg, sign, x_e,
                                                    cheb_out_prec(n_e, x_e));
                if (out) return out;
            }
        }
    }
#endif

    /* Exact integer order with a symbolic / exact-number x: the explicit
     * polynomial (exact for an exact x). Too large stays symbolic. */
    if (n_is_int) {
        if (ni == -1) return expr_new_integer(0);                        /* U_{-1}(x)=0 */
        unsigned long deg; int sign;
        if (!cheb_reflect(ni, &deg, &sign)) return NULL;                 /* absurd order */
        if (deg > CHEB_POLY_CAP) {
            Expr* v = cheb_exact_rational(1, deg, x_e);
            if (v && sign < 0)
                v = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                      (Expr*[]){ expr_new_integer(-1), v }, 2));
            return v;
        }
        return cheb_poly(deg, sign, x_e);
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

Expr* builtin_chebyshev_u(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return builtin_arg_error("ChebyshevU", argc, 2, 2);
    Expr** args = res->data.function.args;
    return chebyshev_u_two_arg(args[0], args[1]);
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

void chebyshev_u_init(void) {
    symtab_add_builtin("ChebyshevU", builtin_chebyshev_u);
    symtab_get_def("ChebyshevU")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

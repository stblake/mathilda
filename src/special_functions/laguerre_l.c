/* Mathilda -- LaguerreL, the (generalized) Laguerre polynomial.
 *
 *   LaguerreL[n, x]      L_n(x)       == LaguerreL[n, 0, x]
 *   LaguerreL[n, a, x]   L_n^(a)(x)   generalized / associated Laguerre
 *
 * The polynomials are orthogonal with weight x^a e^-x and satisfy the ODE
 *   x y'' + (a + 1 - x) y' + n y = 0.
 * They have the closed forms
 *   L_n^(a)(x) = Sum_{i=0}^n (-1)^i Binomial[n+a, n-i] x^i / i!
 *              = Binomial[n+a, n] Hypergeometric1F1[-n, a+1, x]
 * (a = 0 recovers the ordinary L_n).
 *
 * Evaluation is layered so each argument shape takes the cheapest exact or
 * numeric route (mirroring the ChebyshevU / LegendreP design):
 *
 *   exact non-negative integer order n
 *     some inexact argument  ->  the stable three-term recurrence
 *                                  (k+1) L_{k+1} = (2k+1+a-x) L_k - (k+a) L_{k-1}
 *                                evaluated in machine / MPFR complex arithmetic,
 *                                NOT by substituting into the explicit polynomial
 *                                (whose cancellation ruins accuracy for large n,
 *                                 e.g. LaguerreL[100,0,N[5,20]]).
 *     symbolic / exact x, numeric a  ->  the explicit degree-n polynomial in x
 *                                with exact rational coefficients, from the same
 *                                recurrence on coefficient vectors.
 *     symbolic parameter a    ->  the explicit polynomial in a and x, built from
 *                                the sum form and expanded.
 *   non-integer / complex order, some inexact argument
 *                            ->  the confluent-hypergeometric closed form, built
 *                                as an expression over Hypergeometric1F1 so it
 *                                inherits its machine + MPFR, real + complex
 *                                numerics. LaguerreL is entire in x, so a single
 *                                expression covers the whole plane.
 *   special values
 *     L_0(x) = 1,  L_n^(a)(0) = Binomial[n+a, n]  (= 1 when a = 0).
 *   everything else          ->  stays symbolic (return NULL).
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * The machine-precision double / double-complex recurrence lives in
 * sf_machine.c (sf_machine_laguerre_l[_assoc][_c]) and is shared bit-for-bit
 * with the NDArray element-wise kernel and the Compile[] VM (ndkernels.c). The
 * scalar machine path here dispatches to the same helpers.
 *
 * Deferred (left symbolic, matching ChebyshevU): symbolic Series /
 * SeriesCoefficient closed forms, the Series expansion at infinity, D[] of a
 * symbolic order (the D[LaguerreL[...],x] rules live in calculus/deriv.c), and
 * bespoke tight Interval enclosures (the generic NumericFunction + D[] interval
 * threading in eval.c handles those best-effort). Negative integer order with an
 * inexact argument also stays symbolic.
 */
#include "laguerre_l.h"

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "arithmetic.h"      /* make_complex, is_complex */
#include "attr.h"
#include "common.h"          /* builtin_arg_error */
#include "eval.h"            /* eval_and_free */
#include "expr.h"
#include "ndarray.h"         /* is_ndarray, ndt_get/_set, expr_new_ndarray_like */
#include "numeric.h"         /* numeric_min_inexact_bits, get_approx_mpfr */
#include "numeric_complex.h" /* ncpx toolkit */
#include "sf_machine.h"      /* sf_machine_laguerre_l[_assoc][_c] */
#include "symtab.h"
#include "../sym_names.h"

#ifdef USE_MPFR
#include <mpfr.h>
#define PRND MPFR_RNDN
#endif

/* Largest order for which the explicit rational-coefficient polynomial is built
 * (a = 0 / numeric a); beyond this a symbolic/exact-x call stays symbolic. */
#define LAG_POLY_CAP 2000
/* Largest order for which the symbolic-parameter (a is a symbol) polynomial is
 * expanded; the a-coefficients are themselves degree-(n-i) polynomials, so the
 * expansion is far heavier than the numeric case. */
#define LAG_SYM_POLY_CAP 200
/* Largest |order| for which the numeric recurrence runs. */
#define LAG_NUMERIC_DEG_CAP 2000000

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
/* Argument-shape predicates / extractors                             */
/* ------------------------------------------------------------------ */

/* Recognise an exact integer (machine or GMP) argument that fits a long. */
static bool lag_exact_int(const Expr* e, int64_t* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = e->data.integer; return true; }
    if (e->type == EXPR_BIGINT) {
        if (!mpz_fits_slong_p(e->data.bigint)) return false;
        *out = (int64_t)mpz_get_si(e->data.bigint);
        return true;
    }
    return false;
}

/* True if `e` is exactly the integer `v`. */
static bool lag_is_int(const Expr* e, int64_t v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool lag_is_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if `e` carries an inexact numeric leaf (also inside Complex[..]). */
static bool lag_is_inexact(const Expr* e) {
    if (lag_is_inexact_leaf(e)) return true;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im))
        return lag_is_inexact_leaf(re) || lag_is_inexact_leaf(im);
    return false;
}

/* True iff `e` is a pure numeric quantity (Integer/Real/MPFR/BigInt, or
 * Rational/Complex of such) -- i.e. the parameter a is a concrete number. */
static bool lag_is_number(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_REAL:
        case EXPR_BIGINT:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        default: break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        (e->data.function.head->data.symbol.name == SYM_Rational ||
         e->data.function.head->data.symbol.name == SYM_Complex)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!lag_is_number(e->data.function.args[i])) return false;
        return true;
    }
    return false;
}

/* Exact rational/integer value of `e` into `out` (Integer/BigInt/Rational). */
static bool lag_mpq_from_expr(const Expr* e, mpq_t out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) {
        mpq_set_si(out, (long)e->data.integer, 1);
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        mpq_set_z(out, e->data.bigint);
        return true;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational &&
        e->data.function.arg_count == 2) {
        mpz_t num, den;
        mpz_init(num); mpz_init(den);
        expr_to_mpz(e->data.function.args[0], num);   /* Rational args are integer-like */
        expr_to_mpz(e->data.function.args[1], den);
        mpq_set_num(out, num);
        mpq_set_den(out, den);
        mpq_canonicalize(out);
        mpz_clear(num); mpz_clear(den);
        return true;
    }
    return false;
}

/* Machine-double value of a real/integer leaf (declines MPFR/BigInt/complex). */
static bool lag_leaf_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    return false;
}

/* If `x_e` is a machine real, or a Complex of machine real/integer parts, fill
 * (*xr, *xi) and set *real; decline anything carrying MPFR/BigInt.
 * Only the MPFR numeric path below calls this, so its definition is guarded to
 * match -- otherwise it is an unused static under USE_MPFR=0 (-Werror). */
#ifdef USE_MPFR
static bool lag_machine_xy(const Expr* x_e, double* xr, double* xi, bool* real) {
    if (x_e && x_e->type == EXPR_REAL) {
        *xr = x_e->data.real; *xi = 0.0; *real = true; return true;
    }
    if (x_e && x_e->type == EXPR_INTEGER) {
        *xr = (double)x_e->data.integer; *xi = 0.0; *real = true; return true;
    }
    Expr *re, *im;
    if (is_complex((Expr*)x_e, &re, &im) &&
        lag_leaf_double(re, xr) && lag_leaf_double(im, xi)) {
        *real = (*xi == 0.0);
        return true;
    }
    return false;
}
#endif /* USE_MPFR */

#ifdef USE_MPFR
/* True iff the numeric result Expr `e` actually collapsed to a number. */
static bool lag_is_number_leaf(const Expr* e) {
    return lag_is_number(e);
}

/* Output precision: the minimum precision among the inexact arguments, floored
 * at 53 bits. `a_e` may be NULL (the two-argument form). */
static mpfr_prec_t lag_out_prec(const Expr* n_e, const Expr* a_e, const Expr* x_e) {
    long best = 0;
    const Expr* args[3] = { n_e, a_e, x_e };
    for (int i = 0; i < 3; i++) {
        if (args[i] && lag_is_inexact(args[i])) {
            long b = numeric_min_inexact_bits(args[i]);
            if (b > 0) best = (best == 0 || b < best) ? b : best;
        }
    }
    if (best < 53) best = 53;
    return (mpfr_prec_t)best;
}

/* If `n_e` is an inexact real (or MPFR) whose value is an integer within the
 * cap, return true and set *k (so an inexact integer-valued order evaluates via
 * the stable recurrence rather than the general closed form). */
static bool lag_inexact_int_order(const Expr* n_e, int64_t* k) {
    double v;
    if (n_e->type == EXPR_REAL) v = n_e->data.real;
    else if (n_e->type == EXPR_MPFR) v = mpfr_get_d(n_e->data.mpfr, PRND);
    else return false;
    if (v != floor(v)) return false;
    if (v < -(double)LAG_NUMERIC_DEG_CAP || v > (double)LAG_NUMERIC_DEG_CAP)
        return false;
    *k = (int64_t)v;
    return true;
}

/* Build a numeric result Expr from an (re, im) MPFR pair at out_prec bits. */
static Expr* lag_make_result(const mpfr_t re, const mpfr_t im,
                             mpfr_prec_t out_prec, bool force_real) {
    bool real = force_real || mpfr_zero_p(im);
    if (out_prec <= 53) {
        double dr = mpfr_get_d(re, PRND);
        if (real) {
            if (isinf(dr) && !mpfr_inf_p(re)) return expr_new_mpfr_copy(re);
            return expr_new_real(dr);
        }
        double di = mpfr_get_d(im, PRND);
        bool overflow = (isinf(dr) && !mpfr_inf_p(re)) ||
                        (isinf(di) && !mpfr_inf_p(im));
        if (overflow)
            return make_complex(expr_new_mpfr_copy(re), expr_new_mpfr_copy(im));
        return make_complex(expr_new_real(dr), expr_new_real(di));
    }
    Expr* rr = expr_new_mpfr_bits(out_prec);
    mpfr_set(rr->data.mpfr, re, PRND);
    if (real) return rr;
    Expr* ii = expr_new_mpfr_bits(out_prec);
    mpfr_set(ii->data.mpfr, im, PRND);
    return make_complex(rr, ii);
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* mpq coefficient arrays + canonical Expr conversion                 */
/* ------------------------------------------------------------------ */

static mpq_t* lag_mpq_array_new(size_t len) {
    mpq_t* a = (mpq_t*)malloc(len * sizeof(mpq_t));
    if (!a) return NULL;
    for (size_t i = 0; i < len; i++) mpq_init(a[i]);   /* value 0 */
    return a;
}

static void lag_mpq_array_free(mpq_t* a, size_t len) {
    if (!a) return;
    for (size_t i = 0; i < len; i++) mpq_clear(a[i]);
    free(a);
}

static Expr* lag_expr_from_mpz(const mpz_t z) {
    return expr_bigint_normalize(expr_new_bigint_from_mpz(z));
}

/* Build a canonical Integer/BigInt/Rational from an mpq. */
static Expr* lag_expr_from_mpq(const mpq_t q) {
    mpz_t num, den;
    mpz_init(num); mpz_init(den);
    mpq_get_num(num, q);
    mpq_get_den(den, q);
    Expr* en = lag_expr_from_mpz(num);
    Expr* out;
    if (mpz_cmp_ui(den, 1) == 0) {
        out = en;
    } else {
        Expr* ed = lag_expr_from_mpz(den);
        Expr* inv = mk_fn2("Power", ed, expr_new_integer(-1));
        out = eval_and_free(mk_fn2("Times", en, inv));
    }
    mpz_clear(num); mpz_clear(den);
    return out;
}

/* Build Sum_i c[i] x^i (i = 0..deg), skipping zero coefficients. */
static Expr* lag_poly_from_coeffs(mpq_t* c, size_t deg, const Expr* x) {
    Expr** terms = (Expr**)malloc((deg + 1) * sizeof(Expr*));
    if (!terms) return NULL;
    size_t count = 0;
    for (size_t i = 0; i <= deg; i++) {
        if (mpq_sgn(c[i]) == 0) continue;
        Expr* coeff = lag_expr_from_mpq(c[i]);
        Expr* term;
        if (i == 0) {
            term = coeff;
        } else if (i == 1) {
            term = mk_fn2("Times", coeff, expr_copy((Expr*)x));
        } else {
            Expr* xi = mk_fn2("Power", expr_copy((Expr*)x),
                              expr_new_integer((int64_t)i));
            term = mk_fn2("Times", coeff, xi);
        }
        terms[count++] = term;
    }
    Expr* poly;
    if (count == 0)      poly = expr_new_integer(0);
    else if (count == 1) poly = terms[0];
    else                 poly = expr_new_function(expr_new_symbol("Plus"),
                                                  terms, count);
    free(terms);
    return eval_and_free(poly);
}

/* ------------------------------------------------------------------ */
/* Explicit polynomial                                                 */
/* ------------------------------------------------------------------ */

/* L_n^(a)(x) for an exact / absent numeric parameter a, by the three-term
 * recurrence carried on coefficient vectors in x:
 *   L_0 = 1,  L_1 = (1 + a) - x,
 *   L_{k+1}[i] = ((2k+1+a) L_k[i] - L_k[i-1] - (k+a) L_{k-1}[i]) / (k+1).
 * Exact rational arithmetic, so no cancellation. */
static Expr* lag_poly_numeric_a(unsigned long n, const mpq_t a, const Expr* x) {
    /* Closed form, O(n) coefficient operations instead of the recurrence's
     * O(n^2):  c_i = (-1)^i binom(n + a, n - i) / i!, so
     *     c_0 = prod_{j=1..n} (a + j) / j,
     *     c_{i+1} = -c_i (n - i) / ((i + 1)(a + i + 1)).
     * The ratio divides by a + i + 1, so it needs a + j != 0 for j = 1..n;
     * a negative integer in [-n, -1] takes the recurrence below.  (For n = 1000
     * the recurrence's 10^6 operations on ~3000-digit rationals took seconds,
     * and D[LaguerreL[1000, x], x] -- FindRoot's derivative -- paid for it.) */
    bool closed = true;
    if (mpz_cmp_ui(mpq_denref(a), 1) == 0 && mpz_sgn(mpq_numref(a)) < 0
        && mpz_cmpabs_ui(mpq_numref(a), n) <= 0)
        closed = false;
    if (closed) {
        mpq_t* c = lag_mpq_array_new(n + 1);
        if (!c) return NULL;
        mpq_t t;
        mpq_init(t);
        mpq_set_ui(c[0], 1, 1);
        for (unsigned long j = 1; j <= n; j++) {          /* c_0 */
            mpq_set_ui(t, j, 1);
            mpq_add(t, t, a);                             /* a + j        */
            mpq_mul(c[0], c[0], t);
            mpq_set_ui(t, j, 1);
            mpq_div(c[0], c[0], t);
        }
        for (unsigned long i = 0; i < n; i++) {
            mpq_set_ui(t, i + 1, 1);
            mpq_add(t, t, a);                             /* a + i + 1    */
            mpz_mul_ui(mpq_numref(t), mpq_numref(t), i + 1);
            mpq_canonicalize(t);                          /* (i+1)(a+i+1) */
            mpq_div(c[i + 1], c[i], t);
            mpz_mul_ui(mpq_numref(c[i + 1]), mpq_numref(c[i + 1]), n - i);
            mpq_canonicalize(c[i + 1]);
            mpq_neg(c[i + 1], c[i + 1]);
        }
        mpq_clear(t);
        Expr* poly = lag_poly_from_coeffs(c, n, x);
        lag_mpq_array_free(c, n + 1);
        return poly;
    }

    mpq_t* prev = lag_mpq_array_new(n + 1);   /* L_{k-1} */
    mpq_t* cur  = lag_mpq_array_new(n + 1);   /* L_k     */
    mpq_t* nw   = lag_mpq_array_new(n + 1);   /* L_{k+1} */
    if (!prev || !cur || !nw) {
        lag_mpq_array_free(prev, n + 1);
        lag_mpq_array_free(cur, n + 1);
        lag_mpq_array_free(nw, n + 1);
        return NULL;
    }

    mpq_set_ui(prev[0], 1, 1);                 /* L_0 = 1 */
    Expr* poly;
    if (n == 0) {
        poly = lag_poly_from_coeffs(prev, n, x);
    } else {
        mpq_set_ui(cur[0], 1, 1); mpq_add(cur[0], cur[0], a);   /* 1 + a */
        mpq_set_si(cur[1], -1, 1);                               /* -x    */
        if (n == 1) {
            poly = lag_poly_from_coeffs(cur, n, x);
        } else {
            mpq_t c1, c2, t, tmp, kq;
            mpq_inits(c1, c2, t, tmp, kq, (mpq_ptr)0);
            for (unsigned long k = 1; k < n; k++) {
                mpq_set_ui(c1, 2 * k + 1, 1); mpq_add(c1, c1, a);  /* 2k+1+a */
                mpq_set_ui(c2, k, 1);         mpq_add(c2, c2, a);  /* k+a    */
                mpq_set_ui(kq, k + 1, 1);
                for (unsigned long i = 0; i <= n; i++) {
                    mpq_set_ui(t, 0, 1);
                    mpq_mul(tmp, c1, cur[i]); mpq_add(t, t, tmp);     /* (2k+1+a) L_k[i]   */
                    if (i >= 1) mpq_sub(t, t, cur[i - 1]);            /* - L_k[i-1]        */
                    mpq_mul(tmp, c2, prev[i]); mpq_sub(t, t, tmp);    /* - (k+a) L_{k-1}[i]*/
                    mpq_div(t, t, kq);                                /* / (k+1)           */
                    mpq_set(nw[i], t);
                }
                for (unsigned long i = 0; i <= n; i++) {
                    mpq_set(prev[i], cur[i]);
                    mpq_set(cur[i], nw[i]);
                }
            }
            mpq_clears(c1, c2, t, tmp, kq, (mpq_ptr)0);
            poly = lag_poly_from_coeffs(cur, n, x);
        }
    }
    lag_mpq_array_free(prev, n + 1);
    lag_mpq_array_free(cur, n + 1);
    lag_mpq_array_free(nw, n + 1);
    return poly;
}

/* L_n^(a)(x) for a symbolic parameter a, from the sum form
 *   Sum_{i=0}^n (-1)^i / (i! (n-i)!) * Prod_{j=i+1}^n (a + j) * x^i,
 * built as an expression and expanded so it prints as a polynomial in a and x.
 * (x may be an exact or inexact number, giving a numeric-coefficient polynomial
 * in a.) */
static Expr* lag_poly_symbolic_a(unsigned long n, const Expr* a_e, const Expr* x) {
    Expr** terms = (Expr**)malloc((n + 1) * sizeof(Expr*));
    if (!terms) return NULL;
    size_t count = 0;
    mpz_t fi, fni, den;
    mpz_inits(fi, fni, den, (mpz_ptr)0);
    for (unsigned long i = 0; i <= n; i++) {
        mpz_fac_ui(fi, i);          /* i!     */
        mpz_fac_ui(fni, n - i);     /* (n-i)! */
        mpz_mul(den, fi, fni);      /* i! (n-i)! */
        mpq_t cq; mpq_init(cq);
        mpz_set_si(mpq_numref(cq), (i % 2 == 0) ? 1 : -1);
        mpz_set(mpq_denref(cq), den);
        mpq_canonicalize(cq);
        Expr* coeff = lag_expr_from_mpq(cq);
        mpq_clear(cq);

        Expr* prod = NULL;                      /* Prod_{j=i+1}^n (a + j) */
        for (unsigned long j = i + 1; j <= n; j++) {
            Expr* fac = mk_fn2("Plus", expr_copy((Expr*)a_e),
                               expr_new_integer((int64_t)j));
            prod = prod ? mk_fn2("Times", prod, fac) : fac;
        }

        Expr* xp = NULL;                        /* x^i */
        if (i == 1)      xp = expr_copy((Expr*)x);
        else if (i >= 2) xp = mk_fn2("Power", expr_copy((Expr*)x),
                                     expr_new_integer((int64_t)i));

        Expr* term = coeff;
        if (prod) term = mk_fn2("Times", term, prod);
        if (xp)   term = mk_fn2("Times", term, xp);
        terms[count++] = term;
    }
    mpz_clears(fi, fni, den, (mpz_ptr)0);

    Expr* sum = (count == 1) ? terms[0]
                             : expr_new_function(expr_new_symbol("Plus"),
                                                 terms, count);
    free(terms);
    return eval_and_free(mk_fn1("Expand", sum));
}

/* ------------------------------------------------------------------ */
/* Numeric evaluation                                                  */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* L_n^(a)(x) for an inexact argument / parameter and a non-negative integer
 * order n, by the stable forward recurrence in complex MPFR arithmetic. `a_e`
 * may be NULL (a = 0). Returns NULL if an argument is not numeric. */
static Expr* lag_numeric_recurrence(int64_t n, const Expr* a_e, const Expr* x_e,
                                    mpfr_prec_t out_prec) {
    if (n < 0 || n > (int64_t)LAG_NUMERIC_DEG_CAP) return NULL;
    mpfr_prec_t wp = (out_prec < 64 ? 64 : out_prec) + 32;

    ncpx a, x, l0, l1, u, t;
    ncpx_init(&a, wp); ncpx_init(&x, wp); ncpx_init(&l0, wp);
    ncpx_init(&l1, wp); ncpx_init(&u, wp); ncpx_init(&t, wp);

    Expr* out = NULL;
    bool inexact = false;
    bool a_ok = true;
    if (a_e) a_ok = get_approx_mpfr(a_e, a.re, a.im, &inexact);
    else     ncpx_set_ui(&a, 0);

    if (a_ok && get_approx_mpfr(x_e, x.re, x.im, &inexact)) {
        bool real_input = mpfr_zero_p(a.im) && mpfr_zero_p(x.im);

        ncpx_set_ui(&l0, 1);              /* L_0 = 1        */
        ncpx_set_ui(&t, 1);
        ncpx_add(&l1, &t, &a);            /* 1 + a          */
        ncpx_sub(&l1, &l1, &x);           /* L_1 = 1 + a - x */

        ncpx* res;
        if      (n == 0) res = &l0;
        else if (n == 1) res = &l1;
        else if (real_input) {
            /* Real a and x (the common case: a root search, a table): the same
             * recurrence on the real parts alone -- one multiply per product
             * instead of a complex product's four, and the state rotated by
             * swaps instead of copies.  ~4x cheaper per step. */
            mpfr_ptr p0 = l0.re, p1 = l1.re, pu = u.re, pt = t.re;
            for (int64_t k = 1; k < n; k++) {
                mpfr_add_ui(pt, a.re, (unsigned long)(2 * k + 1), PRND);
                mpfr_sub(pt, pt, x.re, PRND);                /* 2k+1+a-x   */
                mpfr_mul(pu, pt, p1, PRND);                  /* (..) L_k   */
                mpfr_add_ui(pt, a.re, (unsigned long)k, PRND);
                mpfr_mul(pt, pt, p0, PRND);                  /* (k+a) L_k-1 */
                mpfr_sub(pu, pu, pt, PRND);
                mpfr_div_ui(pu, pu, (unsigned long)(k + 1), PRND);
                mpfr_ptr tmp = p0; p0 = p1; p1 = pu; pu = tmp;  /* rotate  */
            }
            /* p1 holds L_n; move it into l1.re (the three pointers are a
             * permutation of l0.re, l1.re, u.re). */
            if (p1 != l1.re) mpfr_set(l1.re, p1, PRND);
            res = &l1;
        } else {
            for (int64_t k = 1; k < n; k++) {
                /* u = ((2k+1) + a - x) L_k */
                ncpx_set_ui(&t, (unsigned long)(2 * k + 1));
                ncpx_add(&t, &t, &a);
                ncpx_sub(&t, &t, &x);
                ncpx_mul(&u, &t, &l1, wp);
                /* u -= (k + a) L_{k-1} */
                ncpx_set_ui(&t, (unsigned long)k);
                ncpx_add(&t, &t, &a);
                ncpx_mul(&t, &t, &l0, wp);
                ncpx_sub(&u, &u, &t);
                /* u /= (k + 1) */
                mpfr_div_ui(u.re, u.re, (unsigned long)(k + 1), PRND);
                mpfr_div_ui(u.im, u.im, (unsigned long)(k + 1), PRND);
                ncpx_set(&l0, &l1);
                ncpx_set(&l1, &u);
            }
            res = &l1;
        }
        if (!mpfr_nan_p(res->re) && !mpfr_nan_p(res->im))
            out = lag_make_result(res->re, res->im, out_prec, real_input);
    }

    ncpx_clear(&a); ncpx_clear(&x); ncpx_clear(&l0);
    ncpx_clear(&l1); ncpx_clear(&u); ncpx_clear(&t);
    return out;
}
#endif /* USE_MPFR */

/* Non-integer / complex order numeric via the confluent-hypergeometric form,
 * built as an expression so it inherits Hypergeometric1F1's machine / MPFR,
 * real / complex numerics. `a_e` may be NULL (a = 0):
 *   LaguerreL[n, x]    = Hypergeometric1F1[-n, 1, x]
 *   LaguerreL[n, a, x] = Binomial[n+a, n] Hypergeometric1F1[-n, a+1, x].
 * Returns NULL if the result did not collapse to a number. */
static Expr* lag_numeric_hyp(const Expr* n_e, const Expr* a_e, const Expr* x_e) {
    Expr* neg_n = mk_fn2("Times", expr_new_integer(-1), expr_copy((Expr*)n_e));
    Expr* b, *prefactor;
    if (a_e) {
        b = mk_fn2("Plus", expr_copy((Expr*)a_e), expr_new_integer(1));   /* a + 1 */
        Expr* top = mk_fn2("Plus", expr_copy((Expr*)n_e), expr_copy((Expr*)a_e));
        prefactor = mk_fn2("Binomial", top, expr_copy((Expr*)n_e));       /* C[n+a, n] */
    } else {
        b = expr_new_integer(1);
        prefactor = NULL;
    }
    Expr* hyp = mk_fn3("Hypergeometric1F1", neg_n, b, expr_copy((Expr*)x_e));
    Expr* expr = prefactor ? mk_fn2("Times", prefactor, hyp) : hyp;
    Expr* out = eval_and_free(expr);
#ifdef USE_MPFR
    if (out && lag_is_number_leaf(out)) return out;
#else
    if (out && lag_is_number(out)) return out;
#endif
    expr_free(out);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Driver                                                              */
/* ------------------------------------------------------------------ */

/* `a_e` is NULL for the two-argument form LaguerreL[n, x] (a = 0). */
static Expr* laguerre_l_eval(Expr* n_e, Expr* a_e, Expr* x_e) {
    int64_t ni;
    bool n_is_int = lag_exact_int(n_e, &ni);
    bool x_inexact = lag_is_inexact(x_e);
    bool a_inexact = a_e && lag_is_inexact(a_e);
    bool a_numeric = (a_e == NULL) || lag_is_number(a_e);

    /* Special values. */
    if (n_is_int && ni == 0) return expr_new_integer(1);                /* L_0 = 1 */
    if (!n_is_int && lag_is_int(x_e, 0)) {                              /* L_n^(a)(0) */
        if (a_e == NULL || lag_is_int(a_e, 0)) return expr_new_integer(1);
        Expr* top = mk_fn2("Plus", expr_copy(n_e), expr_copy(a_e));
        return eval_and_free(mk_fn2("Binomial", top, expr_copy(n_e)));  /* Binomial[n+a, n] */
    }

#ifdef USE_MPFR
    /* Numeric recurrence: non-negative integer order (exact or inexact-valued),
     * numeric parameter, and some inexact argument. Avoids the catastrophic
     * cancellation of substituting into the explicit polynomial. */
    {
        int64_t k = 0;
        bool have_int_order = n_is_int ? (k = ni, ni >= 0)
                                       : (lag_inexact_int_order(n_e, &k) && k >= 0);
        bool want_numeric = x_inexact || a_inexact || lag_is_inexact(n_e);
        if (have_int_order && want_numeric && a_numeric) {
            mpfr_prec_t op = lag_out_prec(n_e, a_e, x_e);
            /* Machine fast path: the raw-double / double-complex recurrence,
             * bit-identical to the ND kernel and Compile[]. Taken only at
             * machine precision and degree, with a machine-real parameter and a
             * machine real/complex argument; otherwise the MPFR recurrence. */
            if (op <= 53 && k <= 1000) {
                double ad = 0.0;
                bool a_machine = (a_e == NULL) || lag_leaf_double(a_e, &ad);
                double xr, xi; bool xreal;
                if (a_machine && lag_machine_xy(x_e, &xr, &xi, &xreal)) {
                    if (xreal) {
                        double v;
                        if (sf_machine_laguerre_l_assoc((double)k, ad, xr, &v))
                            return expr_new_real(v);
                    } else {
                        double rr, ri;
                        if (sf_machine_laguerre_l_assoc_c((double)k, ad, xr, xi, &rr, &ri))
                            return make_complex(expr_new_real(rr), expr_new_real(ri));
                    }
                }
            }
            Expr* out = lag_numeric_recurrence(k, a_e, x_e, op);
            if (out) return out;
        }
    }
#endif

    /* Explicit polynomial for a non-negative integer order. */
    if (n_is_int && ni >= 0) {
        if (a_numeric) {                       /* a absent or exact number */
            if (ni > LAG_POLY_CAP) return NULL;
            mpq_t aq; mpq_init(aq);
            bool ok = (a_e == NULL) ? (mpq_set_ui(aq, 0, 1), true)
                                    : lag_mpq_from_expr(a_e, aq);
            Expr* poly = ok ? lag_poly_numeric_a((unsigned long)ni, aq, x_e) : NULL;
            mpq_clear(aq);
            if (poly) return poly;
            /* a inexact (not an exact rational) -> fall through to symbolic
             * builder, which substitutes the numeric a and x. */
        }
        /* Reached only when a is a symbol or an inexact number (a == NULL and
         * exact-numeric a both return via the numeric-a path above). */
        if (a_e && ni <= LAG_SYM_POLY_CAP)
            return lag_poly_symbolic_a((unsigned long)ni, a_e, x_e);
        return NULL;
    }

    /* Non-integer / complex order: numeric only when an argument is inexact. */
    if (lag_is_inexact(n_e) || x_inexact || a_inexact) {
        Expr* out = lag_numeric_hyp(n_e, a_e, x_e);
        if (out) return out;
    }
    return NULL;                               /* stay symbolic */
}

/* ------------------------------------------------------------------ */
/* NDArray fast path (machine buffer argument x)                       */
/* ------------------------------------------------------------------ */

/* Thread the machine recurrence helper over a visible NDArray / packed-List
 * argument `x_e`, holding the scalar order n and parameter a constant. This is
 * what makes the GENERALIZED three-argument form fast over arrays: the 2-arg
 * form reaches the binary ND kernel (eval.c step 4b), but a 3-arg call has no
 * arity->=3 REPL element-wise path, so without this it fell through to the
 * builtin and evaluated the explicit degree-n polynomial over the buffer --
 * correct but ~20-40x slower, and not bit-identical to the scalar path.
 *
 * Reuses sf_machine_laguerre_l_assoc[_c] -- the SAME helper the scalar machine
 * path and Compile[] call -- so the result is bit-for-bit identical to
 * scalar == List == Compile for every order the helper accepts (n a
 * non-negative integer <= 1000, finite result). `a_e` may be NULL (a = 0).
 *
 * Returns NULL to decline ("not my case"): the caller then runs the ordinary
 * scalar/symbolic driver, which keeps the established behaviour for an exact
 * (int64) buffer, a symbolic/non-machine a, or a non-integer order. A buffer
 * the helper declines mid-loop (order > 1000, overflow) is delisted and
 * re-evaluated so each element takes the accurate MPFR scalar route. */
static Expr* laguerre_l_ndarray(const Expr* n_e, const Expr* a_e,
                                const Expr* x_e, Expr* res) {
    if (!is_ndarray(x_e)) return NULL;
    NDType dt = x_e->data.ndarray.dtype;
    /* An exact-integer (or boolean) buffer stays on the exact element-wise path
     * for consistency with the 2-arg kernel and the Sin/ChebyshevU siblings. */
    if (dt == NDT_INT64 || dt == NDT_BOOL) return NULL;
    bool cplx = ndt_is_complex(dt);

    int64_t ni;
    if (!lag_exact_int(n_e, &ni) || ni < 0) return NULL;   /* non-int/neg order */
    double ad = 0.0;
    if (a_e != NULL && !lag_leaf_double(a_e, &ad)) return NULL; /* non-machine a */
    /* The helper caps at n <= 1000 (above that the scalar path uses MPFR for
     * accuracy); delist so a large-order buffer takes that same accurate route
     * element by element rather than a silently less accurate machine pass. */
    if (ni > 1000) return ndarray_delist_and_reeval(res);

    size_t sz = ndarray_size(x_e);
    void* out = malloc(ndt_elem_size(dt) * (sz ? sz : 1));
    if (!out) return NULL;
    bool ok = true;
    for (size_t i = 0; i < sz && ok; i++) {
        double xr, xi, rr = 0.0, ri = 0.0;
        ndt_get(x_e->data.ndarray.data, i, dt, &xr, &xi);
        ok = cplx ? sf_machine_laguerre_l_assoc_c((double)ni, ad, xr, xi, &rr, &ri)
                  : sf_machine_laguerre_l_assoc((double)ni, ad, xr, &rr);
        if (ok) ndt_set(out, i, dt, rr, ri);
    }
    if (ok)
        return expr_new_ndarray_like(x_e, x_e->data.ndarray.rank,
                                     x_e->data.ndarray.dims, out, dt);
    /* A declined element (order > 1000, overflow) must not reach the scalar
     * driver below, which would re-expand the polynomial over the buffer. */
    free(out);
    return ndarray_delist_and_reeval(res);
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                 */
/* ------------------------------------------------------------------ */

Expr* builtin_laguerre_l(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    if (argc != 2 && argc != 3)
        return builtin_arg_error("LaguerreL", argc, 2, 3);
    Expr* a_e = (argc == 3) ? args[1] : NULL;
    Expr* x_e = (argc == 3) ? args[2] : args[1];
    Expr* nd = laguerre_l_ndarray(args[0], a_e, x_e, res);
    if (nd) return nd;
    return laguerre_l_eval(args[0], a_e, x_e);
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

void laguerre_l_init(void) {
    symtab_add_builtin("LaguerreL", builtin_laguerre_l);
    symtab_get_def("LaguerreL")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

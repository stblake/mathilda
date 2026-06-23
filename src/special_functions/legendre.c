/* Mathilda -- Legendre polynomials and associated Legendre functions.
 *
 *   LegendreP[n, x]        Legendre polynomial / function P_n(x).
 *   LegendreP[n, m, x]     associated Legendre function P_n^m(x) (type 1).
 *   LegendreP[n, m, a, x]  Legendre function of type a (a in {1, 2, 3}).
 *
 * Evaluation is layered so each argument shape takes the cheapest exact or
 * numeric route:
 *
 *   LegendreP[n, x]
 *     exact integer n            ->  the explicit degree-|n'| polynomial in x
 *                                    (n' = n, or -1-n for n < 0, since
 *                                    P_{-1-n} = P_n) with exact rational
 *                                    coefficients, built from the three-term
 *                                    recurrence; an inexact x then evaluates
 *                                    the monomials numerically.
 *     x == 1                     ->  1 (for any order n).
 *     non-integer n, some arg    ->  numeric Gauss series
 *        inexact                      P_n(x) = 2F1(-n, n+1; 1; (1-x)/2),
 *                                    real or complex, machine or MPFR
 *                                    precision (requires |(1-x)/2| < 1).
 *     everything else            ->  stays symbolic (return NULL).
 *
 *   LegendreP[n, m, x] / [n, m, a, x]   (integer n, integer m >= 0)
 *     type 1 (default, a == 1)   ->  (-1)^m (1-x^2)^(m/2) d^m/dx^m P_n(x)
 *                                    (the Rodrigues derivative form; 0 when
 *                                    m > |n'|).
 *     types 2, 3                 ->  C(x) * R_a(x), where
 *                                    C(x) = 2F1Reg(-n, n+1, 1-m, (1-x)/2)
 *                                    is the (terminating, exact) regularized
 *                                    Gauss polynomial and the prefactor is
 *                                      R_2 = (1+x)^(m/2) (1-x)^(-m/2),
 *                                      R_3 = (1+x)^(m/2) (-1+x)^(-m/2).
 *     non-integer / negative m   ->  stays symbolic (return NULL).
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Deferred (left symbolic): symbolic Series / SeriesCoefficient, D[] rules,
 * the non-integer associated and Legendre-function forms, and analytic
 * continuation of the numeric series for |(1-x)/2| >= 1.
 */
#include "legendre.h"
#include "sym_names.h"

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "arithmetic.h"      /* make_rational, is_complex */
#include "attr.h"
#include "eval.h"            /* eval_and_free */
#include "expr.h"
#include "numeric.h"         /* numeric_min_inexact_bits, get_approx_mpfr */
#include "numeric_complex.h" /* ncpx toolkit */
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#define PRND MPFR_RNDN
#endif

/* Largest order for which the explicit polynomial is built; beyond this the
 * call stays symbolic. The recurrence is O(n^2) in growing big rationals. */
#define LEG_POLY_CAP 2000

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                  */
/* ------------------------------------------------------------------ */

/* Recognise an exact integer (machine or GMP) first/second argument that fits
 * a long. Returns its value in *out. */
static bool leg_exact_int(const Expr* e, int64_t* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = e->data.integer; return true; }
    if (e->type == EXPR_BIGINT) {
        if (!mpz_fits_slong_p(e->data.bigint)) return false;
        *out = (int64_t)mpz_get_si(e->data.bigint);
        return true;
    }
    return false;
}

/* True if `e` is exactly the integer 1. */
static bool leg_is_int_one(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 1;
}

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool leg_is_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if `e` carries an inexact numeric leaf (also looking inside Complex). */
static bool leg_is_inexact(const Expr* e) {
    if (leg_is_inexact_leaf(e)) return true;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im))
        return leg_is_inexact_leaf(re) || leg_is_inexact_leaf(im);
    return false;
}

/* Build a canonical Integer/BigInt from an mpz. */
static Expr* leg_expr_from_mpz(const mpz_t z) {
    return expr_bigint_normalize(expr_new_bigint_from_mpz(z));
}

/* Build a canonical Integer/BigInt/Rational from an mpq via Times[num,
 * Power[den, -1]] so the evaluator normalises arbitrary-size components. */
static Expr* leg_expr_from_mpq(const mpq_t q) {
    mpz_t num, den;
    mpz_init(num); mpz_init(den);
    mpq_get_num(num, q);
    mpq_get_den(den, q);
    Expr* en = leg_expr_from_mpz(num);
    Expr* out;
    if (mpz_cmp_ui(den, 1) == 0) {
        out = en;
    } else {
        Expr* ed = leg_expr_from_mpz(den);
        Expr* inv = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ ed, expr_new_integer(-1) }, 2);
        out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ en, inv }, 2));
    }
    mpz_clear(num); mpz_clear(den);
    return out;
}

/* ------------------------------------------------------------------ */
/* mpq coefficient arrays                                              */
/* ------------------------------------------------------------------ */

static mpq_t* leg_mpq_array_new(size_t len) {
    mpq_t* a = (mpq_t*)malloc(len * sizeof(mpq_t));
    if (!a) return NULL;
    for (size_t i = 0; i < len; i++) { mpq_init(a[i]); }
    return a;
}

static void leg_mpq_array_free(mpq_t* a, size_t len) {
    if (!a) return;
    for (size_t i = 0; i < len; i++) mpq_clear(a[i]);
    free(a);
}

/* Fill `out` (length n+1, mpq_init'd, value 0) with the coefficients of the
 * Legendre polynomial P_n(x) = Sum_i out[i] x^i, via the three-term
 * recurrence  k P_k = (2k-1) x P_{k-1} - (k-1) P_{k-2}. */
static void legendre_p_coeffs(unsigned long n, mpq_t* out) {
    if (n == 0) { mpq_set_ui(out[0], 1, 1); return; }

    mpq_t* prev = leg_mpq_array_new(n + 1);   /* P_{k-2} */
    mpq_t* cur  = leg_mpq_array_new(n + 1);   /* P_{k-1} */
    mpq_t* nw   = leg_mpq_array_new(n + 1);   /* P_k                       */
    if (!prev || !cur || !nw) {
        leg_mpq_array_free(prev, n + 1);
        leg_mpq_array_free(cur, n + 1);
        leg_mpq_array_free(nw, n + 1);
        return;
    }

    mpq_set_ui(prev[0], 1, 1);                /* P_0 = 1 */
    mpq_set_ui(cur[1], 1, 1);                 /* P_1 = x */

    if (n == 1) {
        for (size_t i = 0; i <= n; i++) mpq_set(out[i], cur[i]);
    } else {
        mpq_t a, b, t;
        mpq_inits(a, b, t, (mpq_ptr)0);
        for (unsigned long k = 2; k <= n; k++) {
            /* nw[i] = ((2k-1) cur[i-1] - (k-1) prev[i]) / k */
            for (unsigned long i = 0; i <= n; i++) {
                mpq_set_ui(t, 0, 1);
                if (i >= 1) {
                    mpq_set_ui(a, 2 * k - 1, 1);
                    mpq_mul(a, a, cur[i - 1]);
                    mpq_add(t, t, a);
                }
                mpq_set_ui(b, k - 1, 1);
                mpq_mul(b, b, prev[i]);
                mpq_sub(t, t, b);
                mpq_set_ui(a, k, 1);
                mpq_div(t, t, a);
                mpq_set(nw[i], t);
            }
            /* rotate: prev <- cur, cur <- nw */
            for (unsigned long i = 0; i <= n; i++) {
                mpq_set(prev[i], cur[i]);
                mpq_set(cur[i], nw[i]);
            }
        }
        for (size_t i = 0; i <= n; i++) mpq_set(out[i], cur[i]);
        mpq_clears(a, b, t, (mpq_ptr)0);
    }

    leg_mpq_array_free(prev, n + 1);
    leg_mpq_array_free(cur, n + 1);
    leg_mpq_array_free(nw, n + 1);
}

/* Build the monomial polynomial Sum_i c[i] x^i (i = 0..deg), skipping zero
 * coefficients, and evaluate it once. Returns integer 0 if all coefficients
 * vanish. */
static Expr* leg_poly_from_coeffs(mpq_t* c, size_t deg, const Expr* x) {
    Expr** terms = (Expr**)malloc((deg + 1) * sizeof(Expr*));
    if (!terms) return NULL;
    size_t count = 0;

    for (size_t i = 0; i <= deg; i++) {
        if (mpq_sgn(c[i]) == 0) continue;
        Expr* coeff = leg_expr_from_mpq(c[i]);
        Expr* term;
        if (i == 0) {
            term = coeff;
        } else if (i == 1) {
            term = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ coeff, expr_copy((Expr*)x) }, 2);
        } else {
            Expr* xi = expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ expr_copy((Expr*)x),
                                     expr_new_integer((int64_t)i) }, 2);
            term = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ coeff, xi }, 2);
        }
        terms[count++] = term;
    }

    Expr* poly;
    if (count == 0)      poly = expr_new_integer(0);
    else if (count == 1) poly = terms[0];
    else                 poly = expr_new_function(expr_new_symbol(SYM_Plus),
                                                  terms, count);
    free(terms);
    return eval_and_free(poly);
}

/* ------------------------------------------------------------------ */
/* LegendreP[n, x] -- the polynomial / function                        */
/* ------------------------------------------------------------------ */

/* Map an integer order to its non-negative representative (P_{-1-n} = P_n). */
static unsigned long leg_norm_order(int64_t n) {
    if (n < 0) n = -1 - n;
    return (unsigned long)n;
}

#ifdef USE_MPFR
/* Build a numeric result Expr from an (re, im) MPFR pair at out_prec bits: a
 * real leaf when `force_real` or im rounds to zero, else Complex[..]. */
static Expr* leg_make_result(const mpfr_t re, const mpfr_t im,
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

/* Numeric P_n(x) via the Gauss series 2F1(-n, n+1; 1; (1-x)/2). Returns NULL
 * if either argument is non-numeric or the series does not converge
 * (|(1-x)/2| >= 1). */
static Expr* leg_numeric_p(const Expr* n_e, const Expr* x_e,
                           mpfr_prec_t out_prec) {
    mpfr_prec_t wp = (out_prec < 64 ? 64 : out_prec) + 32;

    ncpx nz, xz, a, b, w, term, sum, t1, t2, t3;
    ncpx_init(&nz, wp); ncpx_init(&xz, wp); ncpx_init(&a, wp);
    ncpx_init(&b, wp); ncpx_init(&w, wp); ncpx_init(&term, wp);
    ncpx_init(&sum, wp); ncpx_init(&t1, wp); ncpx_init(&t2, wp);
    ncpx_init(&t3, wp);

    bool ok = true, inexact = false;
    if (!get_approx_mpfr(n_e, nz.re, nz.im, &inexact)) ok = false;
    if (ok && !get_approx_mpfr(x_e, xz.re, xz.im, &inexact)) ok = false;

    Expr* out = NULL;
    if (ok) {
        bool real_inputs = mpfr_zero_p(nz.im) && mpfr_zero_p(xz.im);
        ncpx_neg(&a, &nz);                       /* a = -n            */
        ncpx_set(&b, &nz);                        /* b = n + 1         */
        mpfr_add_ui(b.re, b.re, 1, PRND);
        ncpx_neg(&w, &xz);                        /* w = 1 - x         */
        mpfr_add_ui(w.re, w.re, 1, PRND);
        mpfr_t half; mpfr_init2(half, wp);
        mpfr_set_d(half, 0.5, PRND);
        ncpx_scale(&w, &w, half);                 /* w = (1 - x)/2     */

        mpfr_t magw; mpfr_init2(magw, wp);
        ncpx_abs(magw, &w);
        bool converges = (mpfr_cmp_ui(magw, 1) < 0);

        if (converges) {
            ncpx_set_ui(&term, 1);
            ncpx_set_ui(&sum, 1);
            mpfr_t eps, at, asum, rec;
            mpfr_init2(eps, wp); mpfr_init2(at, wp);
            mpfr_init2(asum, wp); mpfr_init2(rec, wp);
            mpfr_set_ui(eps, 1, PRND);
            mpfr_div_2si(eps, eps, (long)wp - 6, PRND); /* 2^-(wp-6) */

            long maxit = 4 * (long)wp + 200;
            bool done = false;
            for (long k = 1; k <= maxit; k++) {
                /* term *= (a+k-1)(b+k-1) / (k*k) * w   (c = 1 => c+k-1 = k) */
                ncpx_set(&t1, &a); mpfr_add_si(t1.re, t1.re, k - 1, PRND);
                ncpx_set(&t2, &b); mpfr_add_si(t2.re, t2.re, k - 1, PRND);
                ncpx_mul(&t3, &t1, &t2, wp);
                ncpx_mul(&term, &term, &t3, wp);
                ncpx_mul(&term, &term, &w, wp);
                mpfr_set_si(rec, k, PRND);
                mpfr_mul_si(rec, rec, k, PRND);   /* k*k */
                mpfr_ui_div(rec, 1, rec, PRND);   /* 1/(k*k) */
                ncpx_scale(&term, &term, rec);
                ncpx_add(&sum, &sum, &term);

                ncpx_abs(at, &term);
                ncpx_abs(asum, &sum);
                mpfr_add_ui(asum, asum, 1, PRND);
                mpfr_mul(asum, asum, eps, PRND);
                if (mpfr_cmp(at, asum) <= 0) { done = true; break; }
            }
            mpfr_clears(eps, at, asum, rec, (mpfr_ptr)0);
            if (done && !mpfr_nan_p(sum.re) && !mpfr_nan_p(sum.im))
                out = leg_make_result(sum.re, sum.im, out_prec, real_inputs);
        }
        mpfr_clears(half, magw, (mpfr_ptr)0);
    }

    ncpx_clear(&nz); ncpx_clear(&xz); ncpx_clear(&a); ncpx_clear(&b);
    ncpx_clear(&w); ncpx_clear(&term); ncpx_clear(&sum); ncpx_clear(&t1);
    ncpx_clear(&t2); ncpx_clear(&t3);
    return out;
}

/* Output precision: the minimum precision among the inexact arguments. */
static mpfr_prec_t leg_out_prec(const Expr* n_e, const Expr* x_e) {
    long bn = leg_is_inexact(n_e) ? numeric_min_inexact_bits(n_e) : 0;
    long bx = leg_is_inexact(x_e) ? numeric_min_inexact_bits(x_e) : 0;
    long bits;
    if (bn > 0 && bx > 0) bits = (bn < bx) ? bn : bx;
    else                  bits = (bn > bx) ? bn : bx;
    if (bits < 53) bits = 53;
    return (mpfr_prec_t)bits;
}
#endif /* USE_MPFR */

/* LegendreP[n, x]. */
static Expr* legendre_two_arg(Expr* n_e, Expr* x_e) {
    /* P_n(1) = 1 for any order n. */
    if (leg_is_int_one(x_e)) return expr_new_integer(1);

    int64_t ni;
    if (leg_exact_int(n_e, &ni)) {
        unsigned long n = leg_norm_order(ni);
        if (n > LEG_POLY_CAP) return NULL;        /* too large: stay symbolic */
        mpq_t* c = leg_mpq_array_new(n + 1);
        if (!c) return NULL;
        legendre_p_coeffs(n, c);
        Expr* poly = leg_poly_from_coeffs(c, n, x_e);
        leg_mpq_array_free(c, n + 1);
        return poly;
    }

#ifdef USE_MPFR
    /* Non-integer order: numeric only when an argument is inexact. */
    if (leg_is_inexact(n_e) || leg_is_inexact(x_e)) {
        Expr* out = leg_numeric_p(n_e, x_e, leg_out_prec(n_e, x_e));
        if (out) return out;
    }
#endif
    return NULL;                                  /* stay symbolic */
}

/* ------------------------------------------------------------------ */
/* Associated Legendre functions                                       */
/* ------------------------------------------------------------------ */

/* Type-1 associated Legendre P_n^m(x) for integer n and integer m >= 0:
 *   (-1)^m (1-x^2)^(m/2) d^m/dx^m P_n(x).
 * Returns integer 0 when m > deg(P_n). */
static Expr* legendre_assoc_type1(unsigned long n, unsigned long m,
                                  const Expr* x) {
    if (m > n) return expr_new_integer(0);

    mpq_t* c = leg_mpq_array_new(n + 1);
    if (!c) return NULL;
    legendre_p_coeffs(n, c);

    /* Differentiate m times: dcoeff[j] = c[j+m] * (j+m)!/j!. */
    size_t ddeg = n - m;
    mpq_t* d = leg_mpq_array_new(ddeg + 1);
    if (!d) { leg_mpq_array_free(c, n + 1); return NULL; }
    mpq_t fall; mpq_init(fall);
    for (size_t j = 0; j <= ddeg; j++) {
        /* falling factorial (j+m)(j+m-1)...(j+1) = (j+m)!/j! */
        mpq_set_ui(fall, 1, 1);
        for (size_t r = 1; r <= m; r++) {
            mpq_t f; mpq_init(f);
            mpq_set_ui(f, (unsigned long)(j + r), 1);
            mpq_mul(fall, fall, f);
            mpq_clear(f);
        }
        mpq_mul(d[j], c[j + m], fall);
    }
    mpq_clear(fall);
    leg_mpq_array_free(c, n + 1);

    Expr* poly = leg_poly_from_coeffs(d, ddeg, x);
    leg_mpq_array_free(d, ddeg + 1);
    if (!poly) return NULL;

    /* prefactor = (-1)^m (1-x^2)^(m/2) */
    Expr* one_minus_x2 = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_new_integer(1),
                   expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1),
                                  expr_new_function(expr_new_symbol(SYM_Power),
                                      (Expr*[]){ expr_copy((Expr*)x),
                                                 expr_new_integer(2) }, 2) }, 2)
                 }, 2);
    Expr* pref = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ one_minus_x2, make_rational((int64_t)m, 2) }, 2);

    Expr* sign = expr_new_integer((m & 1u) ? -1 : 1);
    Expr* result = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ sign, pref, poly }, 3);
    return eval_and_free(result);
}

/* The regularized Gauss polynomial C(x) = 2F1Reg(-n, n+1, 1-m, (1-x)/2) for
 * integer n and integer m >= 0. The series terminates at k = n and its
 * regularized 1/Gamma(1-m+k) factor kills every term with k < m, leaving
 *   C(x) = Sum_{k=m}^{n} (-n)_k (n+1)_k / ((k-m)! k!) ((1-x)/2)^k.
 * Returns integer 0 when m > n. */
static Expr* legendre_2f1reg_core(int64_t ord, unsigned long m, const Expr* x) {
    unsigned long n = leg_norm_order(ord);
    if (m > n) return expr_new_integer(0);

    /* w = (1 - x)/2 = Times[1/2, Plus[1, Times[-1, x]]] */
    Expr* w = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ make_rational(1, 2),
                   expr_new_function(expr_new_symbol(SYM_Plus),
                       (Expr*[]){ expr_new_integer(1),
                                  expr_new_function(expr_new_symbol(SYM_Times),
                                      (Expr*[]){ expr_new_integer(-1),
                                                 expr_copy((Expr*)x) }, 2) }, 2)
                 }, 2);

    Expr** terms = (Expr**)malloc((n + 1) * sizeof(Expr*));
    if (!terms) { expr_free(w); return NULL; }
    size_t count = 0;

    mpz_t num, den, fk, fkm, jj;
    mpz_inits(num, den, fk, fkm, jj, (mpz_ptr)0);
    mpq_t coeff; mpq_init(coeff);

    for (unsigned long k = m; k <= n; k++) {
        /* num = (-n)_k (n+1)_k = prod_{j=0}^{k-1} (-n + j)(n + 1 + j) */
        mpz_set_ui(num, 1);
        for (unsigned long j = 0; j < k; j++) {
            mpz_set_si(jj, -(long)n + (long)j);     /* -n + j */
            mpz_mul(num, num, jj);
            mpz_set_si(jj, (long)n + 1 + (long)j);  /* n + 1 + j */
            mpz_mul(num, num, jj);
        }
        /* den = (k-m)! * k! */
        mpz_fac_ui(fkm, (unsigned long)(k - m));
        mpz_fac_ui(fk, (unsigned long)k);
        mpz_mul(den, fkm, fk);

        mpq_set_num(coeff, num);
        mpq_set_den(coeff, den);
        mpq_canonicalize(coeff);
        if (mpq_sgn(coeff) == 0) continue;

        Expr* cexpr = leg_expr_from_mpq(coeff);
        Expr* wk = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy(w),
                                  expr_new_integer((int64_t)k) }, 2);
        terms[count++] = expr_new_function(expr_new_symbol(SYM_Times),
                             (Expr*[]){ cexpr, wk }, 2);
    }

    mpz_clears(num, den, fk, fkm, jj, (mpz_ptr)0);
    mpq_clear(coeff);
    expr_free(w);

    Expr* sum;
    if (count == 0)      sum = expr_new_integer(0);
    else if (count == 1) sum = terms[0];
    else                 sum = expr_new_function(expr_new_symbol(SYM_Plus),
                                                 terms, count);
    free(terms);
    /* Expand so the result is a clean polynomial in x. */
    Expr* expanded = expr_new_function(expr_new_symbol(SYM_Expand),
                         (Expr*[]){ sum }, 1);
    return eval_and_free(expanded);
}

/* Build Power[base_plus, Rational(num, 2)] where base_plus = Plus[c, sign*x]. */
static Expr* leg_pm_power(int c, int sign, int64_t num, const Expr* x) {
    Expr* base = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_new_integer(c),
                   (sign == 1)
                       ? expr_copy((Expr*)x)
                       : expr_new_function(expr_new_symbol(SYM_Times),
                             (Expr*[]){ expr_new_integer(-1),
                                        expr_copy((Expr*)x) }, 2) }, 2);
    return expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ base, make_rational(num, 2) }, 2);
}

/* Types 2 and 3 of the associated Legendre function for integer n, m >= 0. */
static Expr* legendre_assoc_type23(int64_t ord, unsigned long m, int type,
                                   const Expr* x) {
    Expr* core = legendre_2f1reg_core(ord, m, x);
    if (!core) return NULL;

    /* (1+x)^(m/2) */
    Expr* p_plus = leg_pm_power(1, 1, (int64_t)m, x);
    /* (1-x)^(-m/2) for type 2, (-1+x)^(-m/2) for type 3 */
    Expr* p_minus = (type == 2)
        ? leg_pm_power(1, -1, -(int64_t)m, x)    /* (1 - x)^(-m/2) */
        : leg_pm_power(-1, 1, -(int64_t)m, x);   /* (-1 + x)^(-m/2) */

    Expr* result = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ core, p_plus, p_minus }, 3);
    return eval_and_free(result);
}

/* Shared entry for LegendreP[n, m, x] and LegendreP[n, m, a, x]. `type` is the
 * Legendre function type (1 = default). */
static Expr* legendre_assoc(Expr* n_e, Expr* m_e, int type, Expr* x_e) {
    int64_t ni, mi;
    if (!leg_exact_int(n_e, &ni) || !leg_exact_int(m_e, &mi)) return NULL;
    if (mi < 0) return NULL;                       /* m < 0 deferred */
    if (leg_norm_order(ni) > LEG_POLY_CAP) return NULL;
    unsigned long m = (unsigned long)mi;

    if (type == 1)
        return legendre_assoc_type1(leg_norm_order(ni), m, x_e);
    return legendre_assoc_type23(ni, m, type, x_e);
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                 */
/* ------------------------------------------------------------------ */

/* Mathematica-style diagnostic for a wrong argument count; returns NULL so the
 * evaluator leaves the call unevaluated. */
static Expr* legendre_emit_argb(size_t argc) {
    fprintf(stderr,
            "LegendreP::argb: LegendreP called with %zu argument%s; "
            "between 2 and 4 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_legendre_p(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 2) return legendre_two_arg(args[0], args[1]);
    if (argc == 3) return legendre_assoc(args[0], args[1], 1, args[2]);
    if (argc == 4) {
        int64_t a;
        if (!leg_exact_int(args[2], &a)) return NULL;  /* type must be 1/2/3 */
        if (a < 1 || a > 3) return NULL;
        return legendre_assoc(args[0], args[1], (int)a, args[3]);
    }
    return legendre_emit_argb(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

void legendre_init(void) {
    symtab_add_builtin("LegendreP", builtin_legendre_p);
    symtab_get_def("LegendreP")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

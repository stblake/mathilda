/* Mathilda -- JacobiP, the Jacobi polynomial P_n^(a,b)(x).
 *
 *   JacobiP[n, a, b, x]
 *
 * The polynomials are orthogonal on [-1, 1] with weight (1-x)^a (1+x)^b and
 * satisfy (1 - x^2) y'' + (b - a - (a + b + 2) x) y' + n (n + a + b + 1) y = 0.
 * Closed forms (DLMF 18.5.7, 15.9.1):
 *   P_n^(a,b)(x) = Sum_{l=0}^n (n+a+b+1)_l (a+l+1)_{n-l} / (l! (n-l)!) ((x-1)/2)^l
 *                = Binomial[n+a, n] 2F1(-n, n+a+b+1; a+1; (1-x)/2).
 * The first is a polynomial in a and b as well as x (no division by a
 * parameter), so it defines the polynomial for every a, b.
 *
 * Evaluation is layered (mirroring LaguerreL / ChebyshevU):
 *
 *   exact non-negative integer degree n
 *     some inexact argument  ->  the forward three-term recurrence (DLMF 18.9.2)
 *                                at machine precision (sf_machine_jacobi_p, the
 *                                helper the NDArray kernel and Compile[] share)
 *                                or in MPFR (real, or complex via ncpx), NOT the
 *                                explicit polynomial: its coefficients grow like
 *                                4^n and cancel catastrophically near x = 1.
 *     exact a, b, x          ->  the explicit polynomial with exact rational
 *                                coefficients (an exact rational x evaluates it
 *                                directly by Horner's rule in (x-1)/2).
 *     symbolic a or b        ->  the sum form, expanded (degree <= JP_SYM_CAP).
 *   non-integer / complex degree, some inexact argument, or a degenerate
 *   recurrence (a + b a negative integer <= -2)
 *                            ->  Binomial[n+a, n] Hypergeometric2F1[...], built as
 *                                an expression so it inherits 2F1's numerics.
 *   special values           ->  P_0 = 1,  P_n^(a,b)(1) = Binomial[n+a, n].
 *   everything else          ->  stays symbolic (return NULL).
 *
 * a = b = 0 is LegendreP; a = b = -1/2 gives ChebyshevT up to the factor
 * Binomial[n - 1/2, n]; GegenbauerC is the a = b = lambda - 1/2 case.
 *
 * Attributes: Listable, NumericFunction, Protected. The machine helper is
 * registered as a four-argument N-ary kernel (ndkernels.c), which gives the
 * packed / NDArray element-wise path and the Compile[] lowering.
 * D[JacobiP[n, a, b, x], x] = (n + a + b + 1)/2 JacobiP[n-1, a+1, b+1, x] lives
 * in calculus/deriv.c.
 */
#include "jacobi_p.h"

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include <gmp.h>

#include "arithmetic.h"      /* make_complex, is_complex */
#include "attr.h"
#include "common.h"          /* builtin_arg_error */
#include "eval.h"            /* eval_and_free */
#include "expr.h"
#include "numeric.h"         /* numeric_min_inexact_bits, get_approx_mpfr */
#include "numeric_complex.h" /* ncpx toolkit */
#include "sf_machine.h"      /* sf_machine_jacobi_p[_c] */
#include "symtab.h"
#include "../sym_names.h"

#ifdef USE_MPFR
#include <mpfr.h>
#define PRND MPFR_RNDN
#endif

/* Largest degree for which the exact rational-coefficient polynomial is built. */
#define JP_POLY_CAP 2000
/* Largest degree expanded with a symbolic parameter (the coefficients are then
 * polynomials in a and b, and the expansion is much heavier). */
#define JP_SYM_CAP 40
/* Largest degree run through the numeric recurrence. */
#define JP_NUMERIC_DEG_CAP 2000000

/* ------------------------------------------------------------------ */
/* Expression shorthands and argument predicates                       */
/* ------------------------------------------------------------------ */

static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* mk_fn1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}

static bool jp_exact_int(const Expr* e, int64_t* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = e->data.integer; return true; }
    if (e->type == EXPR_BIGINT) {
        if (!mpz_fits_slong_p(e->data.bigint)) return false;
        *out = (int64_t)mpz_get_si(e->data.bigint);
        return true;
    }
    return false;
}

static bool jp_is_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* An inexact numeric leaf, or a Complex[..] with one. */
static bool jp_is_inexact(const Expr* e) {
    if (jp_is_inexact_leaf(e)) return true;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im))
        return jp_is_inexact_leaf(re) || jp_is_inexact_leaf(im);
    return false;
}

/* A concrete number: Integer/Real/MPFR/BigInt, or Rational/Complex of such. */
static bool jp_is_number(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_REAL: case EXPR_BIGINT:
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
            if (!jp_is_number(e->data.function.args[i])) return false;
        return true;
    }
    return false;
}

/* Exact rational value of an Integer / BigInt / Rational into `out`. */
static bool jp_mpq_from_expr(const Expr* e, mpq_t out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { mpq_set_si(out, (long)e->data.integer, 1); return true; }
    if (e->type == EXPR_BIGINT)  { mpq_set_z(out, e->data.bigint); return true; }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational &&
        e->data.function.arg_count == 2) {
        const Expr* p = e->data.function.args[0];
        const Expr* q = e->data.function.args[1];
        if ((p->type != EXPR_INTEGER && p->type != EXPR_BIGINT) ||
            (q->type != EXPR_INTEGER && q->type != EXPR_BIGINT)) return false;
        mpz_t num, den;
        mpz_init(num); mpz_init(den);
        expr_to_mpz((Expr*)p, num);
        expr_to_mpz((Expr*)q, den);
        mpq_set_num(out, num);
        mpq_set_den(out, den);
        mpq_canonicalize(out);
        mpz_clear(num); mpz_clear(den);
        return true;
    }
    return false;
}

#ifdef USE_MPFR
/* Machine-double value of a machine Real / Integer (declines MPFR, bignum). */
static bool jp_leaf_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    return false;
}
#endif

/* ------------------------------------------------------------------ */
/* Exact polynomial                                                    */
/* ------------------------------------------------------------------ */

static Expr* jp_expr_from_mpz(const mpz_t z) {
    return expr_bigint_normalize(expr_new_bigint_from_mpz(z));
}

static Expr* jp_expr_from_mpq(const mpq_t q) {
    if (mpz_cmp_ui(mpq_denref(q), 1) == 0) return jp_expr_from_mpz(mpq_numref(q));
    Expr* num = jp_expr_from_mpz(mpq_numref(q));
    Expr* den = jp_expr_from_mpz(mpq_denref(q));
    return eval_and_free(mk_fn2("Times", num, mk_fn2("Power", den, expr_new_integer(-1))));
}

/* The coefficients d_l (l = 0..n) of P_n^(a,b) in powers of u = (x-1)/2:
 *   d_l = (n+s+1)_l (a+l+1)_{n-l} / (l! (n-l)!),   s = a + b,
 * in O(n) rational operations: the rising factorial (n+s+1)_l from the left,
 * the product R_l = (a+l+1)...(a+n) from the right, and l!(n-l)! = n!/C(n,l). */
static void jp_u_coeffs(unsigned long n, const mpq_t a, const mpq_t b, mpq_t* d) {
    mpq_t t, s1;
    mpq_inits(t, s1, (mpq_ptr)0);
    /* R_l into d[l]: R_n = 1, R_{l-1} = (a + l) R_l. */
    mpq_set_ui(d[n], 1, 1);
    for (unsigned long l = n; l >= 1; l--) {
        mpq_set_ui(t, l, 1);
        mpq_add(t, t, a);
        mpq_mul(d[l - 1], d[l], t);
    }
    /* times (n+s+1)_l / (l! (n-l)!). */
    mpq_add(s1, a, b);
    mpq_set_ui(t, n + 1, 1);
    mpq_add(s1, s1, t);                       /* n + s + 1 */
    mpq_t rise; mpq_init(rise);
    mpq_set_ui(rise, 1, 1);
    mpz_t fl, fnl;
    mpz_inits(fl, fnl, (mpz_ptr)0);
    mpz_fac_ui(fnl, n);                       /* (n - 0)! */
    mpz_set_ui(fl, 1);                        /* 0!       */
    for (unsigned long l = 0; l <= n; l++) {
        if (l > 0) {
            mpq_set_ui(t, l - 1, 1);
            mpq_add(t, t, s1);
            mpq_mul(rise, rise, t);           /* (n+s+1)_l */
            mpz_mul_ui(fl, fl, l);            /* l!        */
            mpz_divexact_ui(fnl, fnl, n - l + 1); /* (n-l)! */
        }
        mpq_mul(d[l], d[l], rise);
        mpq_set_z(t, fl);
        mpz_mul(mpq_numref(t), mpq_numref(t), fnl);
        mpq_div(d[l], d[l], t);
    }
    mpz_clears(fl, fnl, (mpz_ptr)0);
    mpq_clears(t, s1, rise, (mpq_ptr)0);
}

static mpq_t* jp_mpq_array_new(size_t len) {
    mpq_t* a = (mpq_t*)malloc(len * sizeof(mpq_t));
    if (!a) return NULL;
    for (size_t i = 0; i < len; i++) mpq_init(a[i]);
    return a;
}

static void jp_mpq_array_free(mpq_t* a, size_t len) {
    if (!a) return;
    for (size_t i = 0; i < len; i++) mpq_clear(a[i]);
    free(a);
}

/* Build Sum_i c[i] x^i (zero coefficients skipped) and evaluate it. */
static Expr* jp_poly_from_coeffs(mpq_t* c, size_t deg, const Expr* x) {
    Expr** terms = (Expr**)malloc((deg + 1) * sizeof(Expr*));
    if (!terms) return NULL;
    size_t count = 0;
    for (size_t i = 0; i <= deg; i++) {
        if (mpq_sgn(c[i]) == 0) continue;
        Expr* coeff = jp_expr_from_mpq(c[i]);
        if (i == 0) { terms[count++] = coeff; continue; }
        Expr* xp = (i == 1) ? expr_copy((Expr*)x)
                            : mk_fn2("Power", expr_copy((Expr*)x), expr_new_integer((int64_t)i));
        terms[count++] = mk_fn2("Times", coeff, xp);
    }
    Expr* poly;
    if (count == 0)      poly = expr_new_integer(0);
    else if (count == 1) poly = terms[0];
    else                 poly = expr_new_function(expr_new_symbol("Plus"), terms, count);
    free(terms);
    return eval_and_free(poly);
}

/* P_n^(a,b)(x) for exact rational a, b and an exact (non-inexact) x. */
static Expr* jp_poly_exact(unsigned long n, const mpq_t a, const mpq_t b, const Expr* x) {
    mpq_t* d = jp_mpq_array_new(n + 1);
    if (!d) return NULL;
    jp_u_coeffs(n, a, b, d);

    Expr* out = NULL;
    mpq_t xq;
    mpq_init(xq);
    if (jp_mpq_from_expr(x, xq)) {
        /* Exact rational x: Horner in u = (x-1)/2, exactly. */
        mpq_t u, acc, one;
        mpq_inits(u, acc, one, (mpq_ptr)0);
        mpq_set_ui(one, 1, 1);
        mpq_sub(u, xq, one);
        mpq_div_2exp(u, u, 1);
        mpq_set(acc, d[n]);
        for (unsigned long l = n; l >= 1; l--) {
            mpq_mul(acc, acc, u);
            mpq_add(acc, acc, d[l - 1]);
        }
        out = jp_expr_from_mpq(acc);
        mpq_clears(u, acc, one, (mpq_ptr)0);
    } else {
        /* Symbolic x: re-expand Sum d_l ((x-1)/2)^l in powers of x by Horner's
         * rule on coefficient vectors: c <- c (x - 1)/2 + d_l. */
        mpq_t* c = jp_mpq_array_new(n + 1);
        if (c) {
            mpq_t t; mpq_init(t);
            mpq_set(c[0], d[n]);
            for (unsigned long l = n; l >= 1; l--) {
                size_t deg = n - l;               /* current degree of c */
                /* c * (x - 1)/2: new[i] = (c[i-1] - c[i]) / 2. */
                for (size_t i = deg + 1; i >= 1; i--) {
                    mpq_sub(t, c[i - 1], c[i]);
                    mpq_div_2exp(c[i], t, 1);
                }
                mpq_neg(c[0], c[0]);
                mpq_div_2exp(c[0], c[0], 1);
                mpq_add(c[0], c[0], d[l - 1]);
            }
            mpq_clear(t);
            out = jp_poly_from_coeffs(c, n, x);
            jp_mpq_array_free(c, n + 1);
        }
    }
    mpq_clear(xq);
    jp_mpq_array_free(d, n + 1);
    return out;
}

/* P_n^(a,b)(x) with a symbolic (or inexact) parameter: the sum form
 *   Sum_l (n+a+b+1)_l (a+l+1)_{n-l} / (l! (n-l)!) ((x-1)/2)^l
 * with the rising factorials written as products, expanded. */
static Expr* jp_poly_symbolic(unsigned long n, const Expr* a, const Expr* b, const Expr* x) {
    Expr** terms = (Expr**)malloc((n + 1) * sizeof(Expr*));
    if (!terms) return NULL;
    mpz_t fl, fnl;
    mpz_inits(fl, fnl, (mpz_ptr)0);
    for (unsigned long l = 0; l <= n; l++) {
        mpz_fac_ui(fl, l);
        mpz_fac_ui(fnl, n - l);
        mpz_mul(fl, fl, fnl);
        Expr* term = mk_fn2("Power", jp_expr_from_mpz(fl), expr_new_integer(-1));
        for (unsigned long j = 0; j < l; j++)       /* (n+a+b+1)_l */
            term = mk_fn2("Times", term,
                          expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
                              expr_copy((Expr*)a), expr_copy((Expr*)b),
                              expr_new_integer((int64_t)(n + 1 + j)) }, 3));
        for (unsigned long j = l + 1; j <= n; j++)  /* (a+l+1)_{n-l} */
            term = mk_fn2("Times", term,
                          mk_fn2("Plus", expr_copy((Expr*)a), expr_new_integer((int64_t)j)));
        if (l > 0) {
            Expr* u = mk_fn2("Times", mk_fn2("Rational", expr_new_integer(1), expr_new_integer(2)),
                             mk_fn2("Plus", expr_copy((Expr*)x), expr_new_integer(-1)));
            term = mk_fn2("Times", term, mk_fn2("Power", u, expr_new_integer((int64_t)l)));
        }
        terms[l] = term;
    }
    mpz_clears(fl, fnl, (mpz_ptr)0);
    Expr* sum = (n == 0) ? terms[0]
                         : expr_new_function(expr_new_symbol("Plus"), terms, n + 1);
    free(terms);
    return eval_and_free(mk_fn1("Expand", sum));
}

/* ------------------------------------------------------------------ */
/* Numeric evaluation                                                  */
/* ------------------------------------------------------------------ */

/* Binomial[n+a, n] Hypergeometric2F1[-n, n+a+b+1, a+1, (1-x)/2], evaluated;
 * NULL unless it collapses to a number. */
static Expr* jp_numeric_hyp(const Expr* n, const Expr* a, const Expr* b, const Expr* x) {
    Expr* pre = mk_fn2("Binomial", mk_fn2("Plus", expr_copy((Expr*)n), expr_copy((Expr*)a)),
                       expr_copy((Expr*)n));
    Expr* p1 = mk_fn2("Times", expr_new_integer(-1), expr_copy((Expr*)n));
    Expr* p2 = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
                   expr_copy((Expr*)n), expr_copy((Expr*)a), expr_copy((Expr*)b),
                   expr_new_integer(1) }, 4);
    Expr* p3 = mk_fn2("Plus", expr_copy((Expr*)a), expr_new_integer(1));
    Expr* w = mk_fn2("Times", mk_fn2("Rational", expr_new_integer(1), expr_new_integer(2)),
                     mk_fn2("Plus", expr_new_integer(1),
                            mk_fn2("Times", expr_new_integer(-1), expr_copy((Expr*)x))));
    Expr* hyp = expr_new_function(expr_new_symbol("Hypergeometric2F1"),
                                  (Expr*[]){ p1, p2, p3, w }, 4);
    Expr* out = eval_and_free(mk_fn2("Times", pre, hyp));
    if (out && jp_is_number(out)) return out;
    expr_free(out);
    return NULL;
}

#ifdef USE_MPFR
/* Output precision: the minimum precision among the inexact arguments, >= 53. */
static mpfr_prec_t jp_out_prec(Expr* const* args, int n) {
    long best = 0;
    for (int i = 0; i < n; i++) {
        if (!jp_is_inexact(args[i])) continue;
        long bits = numeric_min_inexact_bits(args[i]);
        if (bits > 0 && (best == 0 || bits < best)) best = bits;
    }
    if (best < 53) best = 53;
    return (mpfr_prec_t)best;
}

/* An inexact real whose value is an integer in range: its degree in *k. */
static bool jp_inexact_int_order(const Expr* e, int64_t* k) {
    double v;
    if (e->type == EXPR_REAL) v = e->data.real;
    else if (e->type == EXPR_MPFR) v = mpfr_get_d(e->data.mpfr, PRND);
    else return false;
    if (v != floor(v) || v < 0.0 || v > (double)JP_NUMERIC_DEG_CAP) return false;
    *k = (int64_t)v;
    return true;
}

/* A numeric result at out_prec bits: a machine Real/Complex at <= 53 bits
 * (MPFR kept when the value overflows a double), else MPFR leaves. */
static Expr* jp_make_result(const mpfr_t re, const mpfr_t im, mpfr_prec_t out_prec,
                            bool force_real) {
    bool real = force_real || mpfr_zero_p(im);
    if (out_prec <= 53) {
        double dr = mpfr_get_d(re, PRND);
        if (real) {
            if (isinf(dr) && !mpfr_inf_p(re)) return expr_new_mpfr_copy(re);
            return expr_new_real(dr);
        }
        double di = mpfr_get_d(im, PRND);
        if ((isinf(dr) && !mpfr_inf_p(re)) || (isinf(di) && !mpfr_inf_p(im)))
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

/* Bit length of n (guard bits for the recurrence's rounding growth). */
static mpfr_prec_t jp_bitlen(int64_t n) {
    mpfr_prec_t b = 1;
    while (n > 1) { n >>= 1; b++; }
    return b;
}

/* P_n^(a,b)(x) by the forward recurrence in MPFR (real parts only when every
 * input is real, else complex ncpx arithmetic).  Returns NULL when an argument
 * is not numeric or the recurrence is degenerate (leading coefficient zero:
 * a + b a negative integer <= -2); *degenerate is set in the latter case. */
static Expr* jp_numeric_recurrence(int64_t n, const Expr* a_e, const Expr* b_e,
                                   const Expr* x_e, mpfr_prec_t out_prec,
                                   bool* degenerate) {
    *degenerate = false;
    mpfr_prec_t wp = (out_prec < 64 ? 64 : out_prec) + 32 + jp_bitlen(n);
    ncpx a, b, x, p0, p1, pk, c2, t, u;
    ncpx_init(&a, wp); ncpx_init(&b, wp); ncpx_init(&x, wp);
    ncpx_init(&p0, wp); ncpx_init(&p1, wp); ncpx_init(&pk, wp);
    ncpx_init(&c2, wp); ncpx_init(&t, wp); ncpx_init(&u, wp);
    mpfr_t s, d, c1, c3, tt;
    mpfr_inits2(wp, s, d, c1, c3, tt, (mpfr_ptr)0);

    Expr* out = NULL;
    bool inexact = false;
    if (!get_approx_mpfr(a_e, a.re, a.im, &inexact) ||
        !get_approx_mpfr(b_e, b.re, b.im, &inexact) ||
        !get_approx_mpfr(x_e, x.re, x.im, &inexact))
        goto done;

    bool real = mpfr_zero_p(a.im) && mpfr_zero_p(b.im) && mpfr_zero_p(x.im);
    bool params_real = mpfr_zero_p(a.im) && mpfr_zero_p(b.im);

    /* P_1 = (a + 1) + (s + 2)(x - 1)/2 */
    ncpx_set_ui(&p0, 1);
    ncpx_add(&t, &a, &b);                      /* s          */
    mpfr_add_ui(t.re, t.re, 2, PRND);          /* s + 2      */
    ncpx_set(&u, &x);
    mpfr_sub_ui(u.re, u.re, 1, PRND);          /* x - 1      */
    ncpx_mul(&p1, &t, &u, wp);
    mpfr_div_2ui(p1.re, p1.re, 1, PRND);
    mpfr_div_2ui(p1.im, p1.im, 1, PRND);
    ncpx_add(&p1, &p1, &a);
    mpfr_add_ui(p1.re, p1.re, 1, PRND);

    ncpx* res = (n == 0) ? &p0 : &p1;
    if (n >= 2 && params_real) {
        /* Real parameters: s, a^2 - b^2 and the three coefficients are real;
         * only x (and so the P_k) may be complex. */
        mpfr_add(s, a.re, b.re, PRND);
        mpfr_sqr(d, a.re, PRND);
        mpfr_sqr(tt, b.re, PRND);
        mpfr_sub(d, d, tt, PRND);                       /* a^2 - b^2 */
        mpfr_ptr P0r = p0.re, P1r = p1.re, PKr = pk.re;
        mpfr_ptr P0i = p0.im, P1i = p1.im, PKi = pk.im;
        mpfr_t tk, q, cr, ci, w1, w2;
        mpfr_inits2(wp, tk, q, cr, ci, w1, w2, (mpfr_ptr)0);
        for (int64_t k = 2; k <= n && !*degenerate; k++) {
            mpfr_add_si(tk, s, 2 * k, PRND);              /* t = 2k + s     */
            mpfr_add_si(c1, s, k, PRND);                  /* k + s          */
            mpfr_mul_si(c1, c1, 2 * k, PRND);
            mpfr_sub_ui(tt, tk, 2, PRND);                 /* t - 2          */
            mpfr_mul(c1, c1, tt, PRND);                   /* 2k(k+s)(t-2)   */
            if (mpfr_zero_p(c1)) { *degenerate = true; break; }
            mpfr_mul(q, tk, tt, PRND);                    /* t(t-2)         */
            mpfr_sub_ui(tt, tk, 1, PRND);                 /* t - 1          */
            mpfr_mul(cr, q, x.re, PRND);
            mpfr_add(cr, cr, d, PRND);
            mpfr_mul(cr, cr, tt, PRND);                   /* Re c2          */
            mpfr_add_si(c3, a.re, k - 1, PRND);
            mpfr_add_si(w1, b.re, k - 1, PRND);
            mpfr_mul(c3, c3, w1, PRND);
            mpfr_mul(c3, c3, tk, PRND);
            mpfr_mul_2ui(c3, c3, 1, PRND);                /* 2(k+a-1)(k+b-1)t */
            if (real) {
                mpfr_mul(PKr, cr, P1r, PRND);
                mpfr_mul(w1, c3, P0r, PRND);
                mpfr_sub(PKr, PKr, w1, PRND);
                mpfr_div(PKr, PKr, c1, PRND);
            } else {
                mpfr_mul(ci, q, x.im, PRND);
                mpfr_mul(ci, ci, tt, PRND);               /* Im c2          */
                /* PK = (c2 P1 - c3 P0) / c1 */
                mpfr_mul(w1, cr, P1r, PRND);
                mpfr_mul(w2, ci, P1i, PRND);
                mpfr_sub(w1, w1, w2, PRND);
                mpfr_mul(w2, c3, P0r, PRND);
                mpfr_sub(w1, w1, w2, PRND);
                mpfr_mul(w2, cr, P1i, PRND);
                mpfr_fma(w2, ci, P1r, w2, PRND);
                mpfr_fms(w2, c3, P0i, w2, PRND);          /* c3 P0i - (..)  */
                mpfr_neg(w2, w2, PRND);
                mpfr_div(PKi, w2, c1, PRND);
                mpfr_div(PKr, w1, c1, PRND);
            }
            /* rotate P0 <- P1 <- PK by pointer swaps */
            mpfr_ptr r0 = P0r; P0r = P1r; P1r = PKr; PKr = r0;
            mpfr_ptr i0 = P0i; P0i = P1i; P1i = PKi; PKi = i0;
        }
        mpfr_clears(tk, q, cr, ci, w1, w2, (mpfr_ptr)0);
        if (!*degenerate) {
            if (P1r != p1.re) mpfr_set(p1.re, P1r, PRND);
            if (P1i != p1.im) mpfr_set(p1.im, P1i, PRND);
            if (real) mpfr_set_zero(p1.im, +1);
        }
        res = &p1;
    } else if (n >= 2) {
        /* Complex parameters: the general complex recurrence. */
        ncpx s_c, d_c, c1c, c3c, q_c;
        ncpx_init(&s_c, wp); ncpx_init(&d_c, wp); ncpx_init(&c1c, wp);
        ncpx_init(&c3c, wp); ncpx_init(&q_c, wp);
        ncpx_add(&s_c, &a, &b);
        ncpx_mul(&d_c, &a, &a, wp);
        ncpx_mul(&t, &b, &b, wp);
        ncpx_sub(&d_c, &d_c, &t);                        /* a^2 - b^2 */
        for (int64_t k = 2; k <= n; k++) {
            ncpx tc;                                     /* t = 2k + s */
            ncpx_init(&tc, wp);
            ncpx_set(&tc, &s_c); mpfr_add_si(tc.re, tc.re, 2 * k, PRND);
            ncpx_set(&c1c, &s_c); mpfr_add_si(c1c.re, c1c.re, k, PRND);
            ncpx_set(&u, &tc); mpfr_sub_ui(u.re, u.re, 2, PRND);       /* t - 2 */
            ncpx_mul(&c1c, &c1c, &u, wp);
            mpfr_mul_si(c1c.re, c1c.re, 2 * k, PRND);
            mpfr_mul_si(c1c.im, c1c.im, 2 * k, PRND);
            if (mpfr_zero_p(c1c.re) && mpfr_zero_p(c1c.im)) {
                *degenerate = true; ncpx_clear(&tc); break;
            }
            ncpx_mul(&q_c, &tc, &u, wp);                 /* t(t-2)       */
            ncpx_mul(&c2, &q_c, &x, wp);
            ncpx_add(&c2, &c2, &d_c);
            ncpx_set(&u, &tc); mpfr_sub_ui(u.re, u.re, 1, PRND);       /* t - 1 */
            ncpx_mul(&c2, &c2, &u, wp);
            ncpx_set(&c3c, &a); mpfr_add_si(c3c.re, c3c.re, k - 1, PRND);
            ncpx_set(&u, &b);   mpfr_add_si(u.re, u.re, k - 1, PRND);
            ncpx_mul(&c3c, &c3c, &u, wp);
            ncpx_mul(&c3c, &c3c, &tc, wp);
            mpfr_mul_2ui(c3c.re, c3c.re, 1, PRND);
            mpfr_mul_2ui(c3c.im, c3c.im, 1, PRND);
            ncpx_mul(&pk, &c2, &p1, wp);
            ncpx_mul(&t, &c3c, &p0, wp);
            ncpx_sub(&pk, &pk, &t);
            ncpx_div(&pk, &pk, &c1c, wp);
            ncpx_set(&p0, &p1);
            ncpx_set(&p1, &pk);
            ncpx_clear(&tc);
        }
        ncpx_clear(&s_c); ncpx_clear(&d_c); ncpx_clear(&c1c);
        ncpx_clear(&c3c); ncpx_clear(&q_c);
        res = &p1;
    }
    if (!*degenerate && !mpfr_nan_p(res->re) && !mpfr_nan_p(res->im))
        out = jp_make_result(res->re, res->im, out_prec, real);

done:
    ncpx_clear(&a); ncpx_clear(&b); ncpx_clear(&x);
    ncpx_clear(&p0); ncpx_clear(&p1); ncpx_clear(&pk);
    ncpx_clear(&c2); ncpx_clear(&t); ncpx_clear(&u);
    mpfr_clears(s, d, c1, c3, tt, (mpfr_ptr)0);
    return out;
}

/* Machine x: a machine Real / Integer, or Complex of such. */
static bool jp_machine_xy(const Expr* x_e, double* xr, double* xi, bool* real) {
    if (jp_leaf_double(x_e, xr)) { *xi = 0.0; *real = true; return true; }
    Expr *re, *im;
    if (is_complex((Expr*)x_e, &re, &im) && jp_leaf_double(re, xr) && jp_leaf_double(im, xi)) {
        *real = false;
        return true;
    }
    return false;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Driver                                                              */
/* ------------------------------------------------------------------ */

static Expr* jacobi_p_eval(Expr* n_e, Expr* a_e, Expr* b_e, Expr* x_e) {
    int64_t ni = 0;
    bool n_is_int = jp_exact_int(n_e, &ni);
    bool params_numeric = jp_is_number(a_e) && jp_is_number(b_e);
    bool any_inexact = jp_is_inexact(n_e) || jp_is_inexact(a_e) ||
                       jp_is_inexact(b_e) || jp_is_inexact(x_e);

    if (n_is_int && ni == 0) return expr_new_integer(1);       /* P_0 = 1 */
    if (n_is_int && ni < 0) return NULL;                       /* not a polynomial */
    if (x_e->type == EXPR_INTEGER && x_e->data.integer == 1 && !any_inexact)
        return eval_and_free(mk_fn2("Binomial",                /* P_n(1)  */
                                    mk_fn2("Plus", expr_copy(n_e), expr_copy(a_e)),
                                    expr_copy(n_e)));

#ifdef USE_MPFR
    if (any_inexact && params_numeric && jp_is_number(x_e)) {
        int64_t k = 0;
        bool have_int = n_is_int ? (k = ni, true) : jp_inexact_int_order(n_e, &k);
        if (have_int && k <= JP_NUMERIC_DEG_CAP) {
            Expr* args[4] = { n_e, a_e, b_e, x_e };
            mpfr_prec_t op = jp_out_prec(args, 4);
            /* Machine fast path: the raw-double helper, bit-identical to the
             * NDArray kernel and Compile[]. */
            double ad, bd, xr, xi;
            bool xreal;
            /* Exact rational parameters (JacobiP[300, 3, 1/2, 0.9998]) enter as
             * their nearest doubles, as they do in the kernel and Compile[]. */
            if (op <= 53 && common_machine_real_value(a_e, &ad) &&
                common_machine_real_value(b_e, &bd) &&
                jp_machine_xy(x_e, &xr, &xi, &xreal)) {
                double rr, ri;
                if (xreal) {
                    if (sf_machine_jacobi_p((double)k, ad, bd, xr, &rr))
                        return expr_new_real(rr);
                } else if (sf_machine_jacobi_p_c((double)k, ad, bd, xr, xi, &rr, &ri)) {
                    return make_complex(expr_new_real(rr), expr_new_real(ri));
                }
            }
            bool degenerate = false;
            Expr* out = jp_numeric_recurrence(k, a_e, b_e, x_e, op, &degenerate);
            if (out) return out;
            if (!degenerate) return NULL;
            /* degenerate recurrence: the hypergeometric form below */
        }
    }
#endif

    if (n_is_int && !any_inexact) {
        mpq_t aq, bq;
        mpq_init(aq); mpq_init(bq);
        Expr* out = NULL;
        bool exact_params = jp_mpq_from_expr(a_e, aq) && jp_mpq_from_expr(b_e, bq);
        if (exact_params && ni <= JP_POLY_CAP)
            out = jp_poly_exact((unsigned long)ni, aq, bq, x_e);
        mpq_clear(aq); mpq_clear(bq);
        if (exact_params) return out;
        if (ni <= JP_SYM_CAP) return jp_poly_symbolic((unsigned long)ni, a_e, b_e, x_e);
        return NULL;
    }
    /* An inexact argument with a symbolic one: the polynomial in the symbols. */
    if (n_is_int && !(params_numeric && jp_is_number(x_e)))
        return ni <= JP_SYM_CAP ? jp_poly_symbolic((unsigned long)ni, a_e, b_e, x_e) : NULL;

    /* Non-integer degree (or a degenerate recurrence): numeric only. */
    if (any_inexact && params_numeric && jp_is_number(n_e) && jp_is_number(x_e))
        return jp_numeric_hyp(n_e, a_e, b_e, x_e);
    return NULL;
}

Expr* builtin_jacobi_p(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 4) return builtin_arg_error("JacobiP", argc, 4, 4);
    Expr** args = res->data.function.args;
    return jacobi_p_eval(args[0], args[1], args[2], args[3]);
}

void jacobi_p_init(void) {
    symtab_add_builtin("JacobiP", builtin_jacobi_p);
    symtab_get_def("JacobiP")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("JacobiP",
        "JacobiP[n, a, b, x]\n"
        "The Jacobi polynomial P_n^(a,b)(x), orthogonal on [-1, 1] with weight\n"
        "(1-x)^a (1+x)^b: Sum_l (n+a+b+1)_l (a+l+1)_(n-l)/(l! (n-l)!) ((x-1)/2)^l\n"
        "= Binomial[n+a, n] Hypergeometric2F1[-n, n+a+b+1, a+1, (1-x)/2].\n"
        "An integer n >= 0 with exact a, b gives the explicit polynomial; numeric\n"
        "arguments evaluate by the stable three-term recurrence at machine or\n"
        "arbitrary precision (real or complex), a non-integer n through 2F1.\n"
        "JacobiP[n, a, b, 1] = Binomial[n+a, n]. Fast paths over packed arrays\n"
        "and in Compile[]. D[JacobiP[n,a,b,x],x] = (n+a+b+1)/2 JacobiP[n-1,a+1,b+1,x].");
}

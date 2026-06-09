/* Mathilda -- the polylogarithm PolyLog.
 *
 *   PolyLog[n, z]     Li_n(z) = Sum_{k>=1} z^k / k^n      (|z| < 1; analytic
 *                     continuation elsewhere, branch cut [1, Infinity))
 *   PolyLog[n, p, z]  Nielsen generalized polylogarithm S_{n,p}(z).  Accepted
 *                     for surface compatibility but left symbolic -- there is
 *                     no closed-form / numeric engine for it here.
 *
 * Evaluation is layered so each kind of argument takes the cheapest exact or
 * fastest numeric route, mirroring src/gamma.c and src/zeta.c:
 *
 *   exact special values (any z) ->  closed forms
 *       PolyLog[n, 0]  = 0
 *       PolyLog[1, z]  = -Log[1 - z]
 *       PolyLog[0, z]  = z/(1 - z)
 *       PolyLog[-m, z] = Eulerian-number rational function   (m >= 1)
 *       PolyLog[n, 1]  = Zeta[n]                              (integer n >= 2)
 *       PolyLog[n, -1] = (2^(1-n) - 1) Zeta[n]                (integer n >= 2)
 *       PolyLog[2, 1/2] = Pi^2/12 - Log[2]^2/2
 *       PolyLog[3, 1/2] = Log[2]^3/6 - Pi^2 Log[2]/12 + 7 Zeta[3]/8
 *   numeric (>= 1 inexact operand, all numeric):
 *       real s, real -1 < z < 1   -> direct real MPFR power series (fast path)
 *       |z| <= 1/2                -> direct power series (complex MPFR)
 *       1/2 < |z|, |ln z| < 2 Pi  -> Jonquiere / zeta expansion (DLMF 25.12.11/12)
 *       otherwise                 -> stays symbolic
 *   everything else -> stays symbolic (return NULL)
 *
 * The zeta expansion needs zeta(s - k) and Gamma(1 - s).  For real order these
 * are evaluated directly with MPFR (mpfr_zeta / mpfr_gamma); for the rare
 * complex-order case they are obtained by evaluating the Zeta / Gamma builtins
 * (which carry their own arbitrary-precision complex kernels).
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "polylog.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "arithmetic.h"   /* is_rational, make_rational, is_complex, make_complex */
#include "numeric.h"      /* numeric_min_inexact_bits */
#include "attr.h"
#include "eval.h"          /* eval_and_free */
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#define PRND MPFR_RNDN
#endif

/* M_PI is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Largest |m| for which we expand PolyLog[-m, z] into its Eulerian rational. */
#define POLYLOG_NEGINT_CAP 400L

/* ------------------------------------------------------------------ */
/* Small predicates / coercions                                        */
/* ------------------------------------------------------------------ */

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool pl_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

#ifdef USE_MPFR
/* True if `e` is an inexact numeric leaf, or a Complex[..] with an inexact
 * part -- i.e. its presence makes a PolyLog call numeric. */
static bool pl_inexact_anywhere(Expr* e) {
    if (pl_is_inexact(e)) return true;
    Expr *re, *im;
    if (is_complex(e, &re, &im)) return pl_is_inexact(re) || pl_is_inexact(im);
    return false;
}
#endif

/* Extract a strictly-exact machine integer (Integer or BigInt fitting a long). */
static bool pl_exact_int(const Expr* e, long* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (long)e->data.integer; return true; }
    if (e->type == EXPR_BIGINT && mpz_fits_slong_p(e->data.bigint)) {
        *out = mpz_get_si(e->data.bigint); return true;
    }
    return false;
}

/* Coerce an exact-or-real leaf to a double (Integer/Real/MPFR/BigInt/Rational).
 * Fails for complex / symbolic. */
static bool pl_to_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { *out = mpfr_get_d(e->data.mpfr, PRND); return true; }
#endif
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
    int64_t n, d;
    if (is_rational(e, &n, &d)) { *out = (double)n / (double)d; return true; }
    return false;
}

/* True if the order argument is (exactly, or inexactly within 1e-9 of) an
 * integer with zero imaginary part, returning that integer in *m. */
static bool pl_integer_order(const Expr* n, long* m) {
    long e;
    if (pl_exact_int(n, &e)) { *m = e; return true; }
    double v;
    if ((n->type == EXPR_REAL
#ifdef USE_MPFR
         || n->type == EXPR_MPFR
#endif
        ) && pl_to_double(n, &v)) {
        double r = floor(v + 0.5);
        if (fabs(v - r) < 1e-9 && fabs(r) < 1e15) { *m = (long)r; return true; }
    }
    return false;
}

/* True if `e` is the exact rational 1/2. */
static bool pl_is_half(const Expr* e) {
    int64_t n, d;
    return is_rational(e, &n, &d) && n == 1 && d == 2;
}

/* ------------------------------------------------------------------ */
/* Exact closed forms (symbolic; reused for inexact via evaluation)    */
/* ------------------------------------------------------------------ */

/* PolyLog[1, z] = -Log[1 - z]. */
static Expr* pl_order1(Expr* z) {
    Expr* omz = expr_new_function(expr_new_symbol("Subtract"),
                    (Expr*[]){ expr_new_integer(1), expr_copy(z) }, 2);
    Expr* lg  = expr_new_function(expr_new_symbol("Log"), &omz, 1);
    Expr* neg = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), lg }, 2);
    return eval_and_free(neg);
}

/* PolyLog[0, z] = z/(1 - z). */
static Expr* pl_order0(Expr* z) {
    Expr* omz = expr_new_function(expr_new_symbol("Subtract"),
                    (Expr*[]){ expr_new_integer(1), expr_copy(z) }, 2);
    Expr* inv = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ omz, expr_new_integer(-1) }, 2);
    Expr* out = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(z), inv }, 2);
    return eval_and_free(out);
}

/* PolyLog[-m, z] for m >= 1, the Eulerian-number rational:
 *   Li_{-m}(z) = (1/(1-z)^{m+1}) Sum_{k=0}^{m-1} A(m, k) z^{m-k},
 * with A(m, k) the Eulerian numbers, A(m,k) = (k+1)A(m-1,k) + (m-k)A(m-1,k-1).
 * Builds the symbolic expression and evaluates it (exact z -> exact rational,
 * symbolic z -> rational function, inexact z -> numeric). NULL past the cap. */
static Expr* pl_neg_order(long m, Expr* z) {
    if (m < 1 || m > POLYLOG_NEGINT_CAP) return NULL;

    /* Eulerian row A(m, 0..m-1) via the triangle recurrence (GMP, exact). */
    mpz_t* A = (mpz_t*)malloc((size_t)m * sizeof(mpz_t));
    if (!A) return NULL;
    for (long k = 0; k < m; k++) mpz_init(A[k]);
    mpz_set_ui(A[0], 1);                          /* row m = 1: A(1,0) = 1 */
    for (long r = 2; r <= m; r++) {
        /* Update row r from row r-1 in place, going high index to low. */
        for (long k = r - 1; k >= 0; k--) {
            mpz_t left;  mpz_init(left);
            if (k < r - 1) {                       /* (k+1) A(r-1, k) */
                mpz_mul_ui(left, A[k], (unsigned long)(k + 1));
            }
            if (k >= 1) {                          /* + (r-k) A(r-1, k-1) */
                mpz_addmul_ui(left, A[k - 1], (unsigned long)(r - k));
            }
            mpz_set(A[k], left);
            mpz_clear(left);
        }
    }

    /* numerator = Sum_{k=0}^{m-1} A(m,k) z^{m-k}. */
    Expr** terms = (Expr**)malloc((size_t)m * sizeof(Expr*));
    if (!terms) { for (long k = 0; k < m; k++) mpz_clear(A[k]); free(A); return NULL; }
    for (long k = 0; k < m; k++) {
        Expr* coeff = expr_bigint_normalize(expr_new_bigint_from_mpz(A[k]));
        long pw = m - k;                           /* exponent m-k, >= 1 */
        Expr* zk = expr_new_function(expr_new_symbol("Power"),
                       (Expr*[]){ expr_copy(z), expr_new_integer(pw) }, 2);
        terms[k] = expr_new_function(expr_new_symbol("Times"),
                       (Expr*[]){ coeff, zk }, 2);
        mpz_clear(A[k]);
    }
    free(A);

    Expr* numer = (m == 1) ? terms[0]
                  : expr_new_function(expr_new_symbol("Plus"), terms, (size_t)m);
    free(terms);

    /* denom^{-1} = (1 - z)^{-(m+1)}. */
    Expr* omz = expr_new_function(expr_new_symbol("Subtract"),
                    (Expr*[]){ expr_new_integer(1), expr_copy(z) }, 2);
    Expr* denom = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ omz, expr_new_integer(-(m + 1)) }, 2);

    Expr* out = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ numer, denom }, 2);
    return eval_and_free(out);
}

/* PolyLog[n, 1] = Zeta[n] (integer n >= 2). */
static Expr* pl_at_one(long n) {
    Expr* z = expr_new_function(expr_new_symbol("Zeta"),
                  (Expr*[]){ expr_new_integer(n) }, 1);
    return eval_and_free(z);
}

/* PolyLog[n, -1] = (2^(1-n) - 1) Zeta[n] = -eta(n) (integer n >= 2). */
static Expr* pl_at_minus_one(long n) {
    /* coefficient 2^(1-n) - 1 = Subtract[Power[2, 1-n], 1]. */
    Expr* pw = expr_new_function(expr_new_symbol("Power"),
                   (Expr*[]){ expr_new_integer(2), expr_new_integer(1 - n) }, 2);
    Expr* coeff = expr_new_function(expr_new_symbol("Subtract"),
                   (Expr*[]){ pw, expr_new_integer(1) }, 2);
    Expr* zeta = expr_new_function(expr_new_symbol("Zeta"),
                   (Expr*[]){ expr_new_integer(n) }, 1);
    Expr* out = expr_new_function(expr_new_symbol("Times"),
                   (Expr*[]){ coeff, zeta }, 2);
    return eval_and_free(out);
}

/* PolyLog[2, 1/2] = Pi^2/12 - Log[2]^2/2. */
static Expr* pl_two_half(void) {
    Expr* pi2 = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_symbol("Pi"), expr_new_integer(2) }, 2);
    Expr* t1 = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ make_rational(1, 12), pi2 }, 2);
    Expr* log2 = expr_new_function(expr_new_symbol("Log"),
                    (Expr*[]){ expr_new_integer(2) }, 1);
    Expr* log2sq = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ log2, expr_new_integer(2) }, 2);
    Expr* t2 = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ make_rational(-1, 2), log2sq }, 2);
    Expr* out = expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ t1, t2 }, 2);
    return eval_and_free(out);
}

/* PolyLog[3, 1/2] = Log[2]^3/6 - Pi^2 Log[2]/12 + 7 Zeta[3]/8. */
static Expr* pl_three_half(void) {
    Expr* log2a = expr_new_function(expr_new_symbol("Log"),
                    (Expr*[]){ expr_new_integer(2) }, 1);
    Expr* log2cube = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ log2a, expr_new_integer(3) }, 2);
    Expr* t1 = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ make_rational(1, 6), log2cube }, 2);

    Expr* pi2 = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_symbol("Pi"), expr_new_integer(2) }, 2);
    Expr* log2b = expr_new_function(expr_new_symbol("Log"),
                    (Expr*[]){ expr_new_integer(2) }, 1);
    Expr* t2 = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ make_rational(-1, 12), pi2, log2b }, 3);

    Expr* zeta3 = expr_new_function(expr_new_symbol("Zeta"),
                    (Expr*[]){ expr_new_integer(3) }, 1);
    Expr* t3 = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ make_rational(7, 8), zeta3 }, 2);

    Expr* out = expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ t1, t2, t3 }, 3);
    return eval_and_free(out);
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ */
/* Complex-MPFR toolkit (pairs of mpfr_t; no MPC library available).   */
/* Mirrors the `zcx` helpers in src/zeta.c -- alias-safe, explicit      */
/* working precision `p`.                                               */
/* ------------------------------------------------------------------ */

typedef struct { mpfr_t re, im; } pcx;

static void pcx_init(pcx* z, mpfr_prec_t p) { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void pcx_clear(pcx* z)               { mpfr_clear(z->re);    mpfr_clear(z->im);    }
static void pcx_set(pcx* d, const pcx* s)   { mpfr_set(d->re, s->re, PRND); mpfr_set(d->im, s->im, PRND); }

static void pcx_add(pcx* out, const pcx* a, const pcx* b) {
    mpfr_add(out->re, a->re, b->re, PRND);
    mpfr_add(out->im, a->im, b->im, PRND);
}

static void pcx_abs(mpfr_t mag, const pcx* z) { mpfr_hypot(mag, z->re, z->im, PRND); }

static void pcx_mul(pcx* out, const pcx* a, const pcx* b, mpfr_prec_t p) {
    mpfr_t ac, bd, ad, bc;
    mpfr_inits2(p, ac, bd, ad, bc, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, PRND);
    mpfr_mul(bd, a->im, b->im, PRND);
    mpfr_mul(ad, a->re, b->im, PRND);
    mpfr_mul(bc, a->im, b->re, PRND);
    mpfr_sub(out->re, ac, bd, PRND);
    mpfr_add(out->im, ad, bc, PRND);
    mpfr_clears(ac, bd, ad, bc, (mpfr_ptr)0);
}

/* out = a * r, r a real scalar. */
static void pcx_mul_r(pcx* out, const pcx* a, const mpfr_t r) {
    mpfr_mul(out->re, a->re, r, PRND);
    mpfr_mul(out->im, a->im, r, PRND);
}

/* out = a / r, r a real scalar. */
static void pcx_div_r(pcx* out, const pcx* a, const mpfr_t r) {
    mpfr_div(out->re, a->re, r, PRND);
    mpfr_div(out->im, a->im, r, PRND);
}

static void pcx_exp(pcx* out, const pcx* z, mpfr_prec_t p) {
    mpfr_t ea, c, s;
    mpfr_inits2(p, ea, c, s, (mpfr_ptr)0);
    mpfr_exp(ea, z->re, PRND);
    mpfr_sin_cos(s, c, z->im, PRND);
    mpfr_mul(out->re, ea, c, PRND);
    mpfr_mul(out->im, ea, s, PRND);
    mpfr_clears(ea, c, s, (mpfr_ptr)0);
}

/* out = Log(z), principal branch: ln|z| + i Arg(z). */
static void pcx_log(pcx* out, const pcx* z, mpfr_prec_t p) {
    mpfr_t mag, ang;
    mpfr_inits2(p, mag, ang, (mpfr_ptr)0);
    mpfr_hypot(mag, z->re, z->im, PRND);
    mpfr_atan2(ang, z->im, z->re, PRND);
    mpfr_log(out->re, mag, PRND);
    mpfr_set(out->im, ang, PRND);
    mpfr_clears(mag, ang, (mpfr_ptr)0);
}

/* out = sin(z) = sin(a)cosh(b) + i cos(a)sinh(b). */
static void pcx_sin(pcx* out, const pcx* z, mpfr_prec_t p) {
    mpfr_t sa, ca, chb, shb;
    mpfr_inits2(p, sa, ca, chb, shb, (mpfr_ptr)0);
    mpfr_sin_cos(sa, ca, z->re, PRND);
    mpfr_sinh_cosh(shb, chb, z->im, PRND);
    mpfr_mul(out->re, sa, chb, PRND);
    mpfr_mul(out->im, ca, shb, PRND);
    mpfr_clears(sa, ca, chb, shb, (mpfr_ptr)0);
}

/* out = base^expo (principal branch) = exp(expo * Log(base)). */
static void pcx_pow(pcx* out, const pcx* base, const pcx* expo, mpfr_prec_t p) {
    pcx lg, prod;
    pcx_init(&lg, p); pcx_init(&prod, p);
    pcx_log(&lg, base, p);
    pcx_mul(&prod, expo, &lg, p);
    pcx_exp(out, &prod, p);
    pcx_clear(&lg); pcx_clear(&prod);
}

/* ------------------------------------------------------------------ */
/* MPFR coercion / result builders                                     */
/* ------------------------------------------------------------------ */

/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool pl_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, PRND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          PRND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        PRND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          PRND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, PRND);
        mpfr_div_si(out, out, (long)d, PRND);
        return true;
    }
    return false;
}

/* Result precision (bits) under Mathematica contagion: the minimum precision
 * among the inexact leaves of `a` and `b`, floored at machine (53). A machine
 * `Real` argument therefore forces a machine-precision result, even alongside
 * a high-precision MPFR argument: PolyLog[N[1/2], 1.3429`50] is machine, not
 * 50-digit. See numeric_min_inexact_bits. */
static mpfr_prec_t pl_out_prec(const Expr* a, const Expr* b) {
    long pa = numeric_min_inexact_bits(a);
    long pb = numeric_min_inexact_bits(b);
    long m  = (pa && pb) ? (pa < pb ? pa : pb) : (pa ? pa : pb);
    return m < 53 ? 53 : (mpfr_prec_t)m;
}

/* Build the numeric result: Real / Complex[Real,Real] at machine precision,
 * MPFR parts at arbitrary precision; `real_only` forces a real result. */
static Expr* pl_result(const mpfr_t re, const mpfr_t im,
                       mpfr_prec_t out_prec, bool real_only) {
    if (out_prec <= 53) {
        Expr* rr = expr_new_real(mpfr_get_d(re, PRND));
        if (real_only) return rr;
        return make_complex(rr, expr_new_real(mpfr_get_d(im, PRND)));
    }
    Expr* rr = expr_new_mpfr_bits(out_prec);
    mpfr_set(rr->data.mpfr, re, PRND);
    if (real_only) return rr;
    Expr* ii = expr_new_mpfr_bits(out_prec);
    mpfr_set(ii->data.mpfr, im, PRND);
    return make_complex(rr, ii);
}

/* Decompose a numeric value into borrowed (re, im); plain real -> (value, zero). */
static void pl_decompose(Expr* v, Expr* zero, Expr** re, Expr** im) {
    Expr *r, *i;
    if (is_complex(v, &r, &i)) { *re = r; *im = i; }
    else { *re = v; *im = zero; }
}

/* Fill an already-init'd pcx (precision wp) from an evaluated numeric Expr. */
static bool pl_set_pcx_from_expr(pcx* out, Expr* e) {
    if (!e) return false;
    Expr *re, *im;
    if (is_complex(e, &re, &im))
        return pl_set_mpfr(out->re, re) && pl_set_mpfr(out->im, im);
    if (pl_set_mpfr(out->re, e)) { mpfr_set_ui(out->im, 0, PRND); return true; }
    return false;
}

/* Gamma(w) into out. Real w uses mpfr_gamma; complex w routes through the
 * Gamma builtin (its arbitrary-precision Spouge kernel). */
static bool pl_gamma_cx(pcx* out, const pcx* w) {
    if (mpfr_zero_p(w->im)) {
        mpfr_gamma(out->re, w->re, PRND);
        mpfr_set_ui(out->im, 0, PRND);
        return mpfr_number_p(out->re);
    }
    Expr* arg = make_complex(expr_new_mpfr_copy(w->re), expr_new_mpfr_copy(w->im));
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol("Gamma"), &arg, 1));
    bool ok = pl_set_pcx_from_expr(out, r);
    expr_free(r);
    return ok && mpfr_number_p(out->re) && mpfr_number_p(out->im);
}

/* Raw complex zeta(w) via the Zeta builtin. The builtin's Euler-Maclaurin
 * kernel is only well-conditioned for Re(w) >= 1/2; callers must reflect first
 * for the left half-plane (see pl_zeta_cx). */
static bool pl_zeta_eval_raw(pcx* out, const pcx* w) {
    Expr* arg = make_complex(expr_new_mpfr_copy(w->re), expr_new_mpfr_copy(w->im));
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol("Zeta"), &arg, 1));
    bool ok = pl_set_pcx_from_expr(out, r);
    expr_free(r);
    return ok && mpfr_number_p(out->re) && mpfr_number_p(out->im);
}

/* zeta(w) into out (precision wp). Real w uses mpfr_zeta (reflection handled
 * internally). Complex w with Re(w) >= 1/2 calls the Zeta builtin directly;
 * for Re(w) < 1/2 the functional equation
 *   zeta(w) = 2^w pi^(w-1) sin(pi w/2) Gamma(1-w) zeta(1-w)
 * is applied -- a product with no cancellation, so it stays accurate where the
 * raw Euler-Maclaurin sum at large negative Re(w) would lose all precision.
 * Returns false on a non-finite / non-numeric result. */
static bool pl_zeta_cx(pcx* out, const pcx* w, mpfr_prec_t wp) {
    if (mpfr_zero_p(w->im)) {
        mpfr_set_ui(out->im, 0, PRND);
        /* For a real argument in the left half-plane, mpfr_zeta runs its own
         * reflection internally, which is ~5x slower than evaluating zeta on
         * the convergent right half-plane.  When w is comfortably away from an
         * integer (so neither the sin(pi w/2) zero nor the zeta(1-w) pole is in
         * play -- the polylog expansion feeds half-integer orders here) we apply
         * the functional equation ourselves, reusing the fast positive-side
         * mpfr_zeta(1-w).  Integers and near-integers fall through to the direct
         * call, which already handles the trivial zeros / pole exactly. */
        mpfr_t frac;
        mpfr_init2(frac, wp);
        mpfr_round(frac, w->re);                       /* nearest integer */
        mpfr_sub(frac, w->re, frac, PRND);             /* distance to it */
        mpfr_abs(frac, frac, PRND);
        bool near_int = mpfr_cmp_d(frac, 1e-4) < 0;
        mpfr_clear(frac);
        if (mpfr_cmp_d(w->re, 0.5) < 0 && !near_int) {
            mpfr_t omw, z1, g, sn, fac, pi, ln2, lnpi, t;
            mpfr_inits2(wp, omw, z1, g, sn, fac, pi, ln2, lnpi, t, (mpfr_ptr)0);
            mpfr_const_pi(pi, PRND);
            mpfr_set_ui(ln2, 2, PRND); mpfr_log(ln2, ln2, PRND);
            mpfr_log(lnpi, pi, PRND);
            mpfr_ui_sub(omw, 1, w->re, PRND);          /* 1 - w (> 1/2) */
            mpfr_zeta(z1, omw, PRND);                   /* fast: right half-plane */
            mpfr_gamma(g, omw, PRND);
            mpfr_mul(sn, w->re, pi, PRND);
            mpfr_div_2ui(sn, sn, 1, PRND);
            mpfr_sin(sn, sn, PRND);                     /* sin(pi w / 2) */
            mpfr_mul(fac, w->re, ln2, PRND);            /* w ln2 */
            mpfr_sub_ui(t, w->re, 1, PRND);
            mpfr_mul(t, t, lnpi, PRND);                 /* (w-1) ln pi */
            mpfr_add(fac, fac, t, PRND);
            mpfr_exp(fac, fac, PRND);                   /* 2^w pi^(w-1) */
            mpfr_mul(out->re, fac, sn, PRND);
            mpfr_mul(out->re, out->re, g, PRND);
            mpfr_mul(out->re, out->re, z1, PRND);
            mpfr_clears(omw, z1, g, sn, fac, pi, ln2, lnpi, t, (mpfr_ptr)0);
            return mpfr_number_p(out->re);
        }
        mpfr_zeta(out->re, w->re, PRND);
        return mpfr_number_p(out->re);
    }
    if (mpfr_cmp_d(w->re, 0.5) >= 0)
        return pl_zeta_eval_raw(out, w);

    /* Reflection for the left half-plane. */
    pcx omw, z1, g, sn, fac, t;
    pcx_init(&omw, wp); pcx_init(&z1, wp); pcx_init(&g, wp);
    pcx_init(&sn, wp); pcx_init(&fac, wp); pcx_init(&t, wp);
    mpfr_t pi, ln2, lnpi, h;
    mpfr_inits2(wp, pi, ln2, lnpi, h, (mpfr_ptr)0);
    mpfr_const_pi(pi, PRND);
    mpfr_set_ui(ln2, 2, PRND); mpfr_log(ln2, ln2, PRND);
    mpfr_log(lnpi, pi, PRND);

    mpfr_ui_sub(omw.re, 1, w->re, PRND);          /* 1 - w */
    mpfr_neg(omw.im, w->im, PRND);

    bool ok = pl_zeta_eval_raw(&z1, &omw) && pl_gamma_cx(&g, &omw);
    if (ok) {
        /* sin(pi w / 2). */
        mpfr_div_ui(h, pi, 2, PRND);
        mpfr_mul(t.re, w->re, h, PRND);
        mpfr_mul(t.im, w->im, h, PRND);
        pcx_sin(&sn, &t, wp);
        /* fac = 2^w pi^(w-1) = exp(w ln2 + (w-1) ln pi). */
        mpfr_mul(t.re, w->re, ln2, PRND);
        mpfr_mul(t.im, w->im, ln2, PRND);
        mpfr_set(fac.re, t.re, PRND); mpfr_set(fac.im, t.im, PRND);
        mpfr_sub_ui(t.re, w->re, 1, PRND);        /* (w-1) re */
        mpfr_mul(t.re, t.re, lnpi, PRND);
        mpfr_mul(t.im, w->im, lnpi, PRND);
        mpfr_add(fac.re, fac.re, t.re, PRND);
        mpfr_add(fac.im, fac.im, t.im, PRND);
        pcx_exp(&fac, &fac, wp);
        /* out = fac * sn * g * z1. */
        pcx_mul(out, &fac, &sn, wp);
        pcx_mul(out, out, &g, wp);
        pcx_mul(out, out, &z1, wp);
        ok = mpfr_number_p(out->re) && mpfr_number_p(out->im);
    }
    mpfr_clears(pi, ln2, lnpi, h, (mpfr_ptr)0);
    pcx_clear(&omw); pcx_clear(&z1); pcx_clear(&g);
    pcx_clear(&sn); pcx_clear(&fac); pcx_clear(&t);
    return ok;
}

/* ------------------------------------------------------------------ */
/* Numeric kernels                                                     */
/* ------------------------------------------------------------------ */

/* Direct real power series Li_s(z) = Sum_{k>=1} z^k k^-s for real s and real
 * -1 < z < 1.  Pure-MPFR fast path (no complex arithmetic). */
static Expr* pl_real_series(const Expr* sre, const Expr* zre,
                            mpfr_prec_t out_prec, mpfr_prec_t wp) {
    mpfr_t s, z, S, zk, term, kln, ks, eps, at, as;
    mpfr_inits2(wp, s, z, S, zk, term, kln, ks, eps, at, as, (mpfr_ptr)0);
    pl_set_mpfr(s, sre);
    pl_set_mpfr(z, zre);
    mpfr_set_ui(S, 0, PRND);
    mpfr_set(zk, z, PRND);                        /* z^1 */
    mpfr_set_ui(eps, 1, PRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 8 ? wp - 6 : 1), PRND);

    for (unsigned long k = 1; k < 50000000UL; k++) {
        if (k == 1) {
            mpfr_set(term, zk, PRND);             /* 1^-s = 1 */
        } else {
            mpfr_set_ui(kln, k, PRND);
            mpfr_log(kln, kln, PRND);
            mpfr_mul(ks, s, kln, PRND);
            mpfr_neg(ks, ks, PRND);
            mpfr_exp(ks, ks, PRND);               /* k^-s */
            mpfr_mul(term, zk, ks, PRND);
        }
        mpfr_add(S, S, term, PRND);
        if (k >= 3) {
            mpfr_abs(at, term, PRND);
            mpfr_abs(as, S, PRND);
            mpfr_mul(as, as, eps, PRND);
            if (mpfr_cmp(at, as) < 0) break;
        }
        mpfr_mul(zk, zk, z, PRND);
    }

    Expr* out;
    if (out_prec <= 53) {
        out = expr_new_real(mpfr_get_d(S, PRND));
    } else {
        out = expr_new_mpfr_bits(out_prec);
        mpfr_set(out->data.mpfr, S, PRND);
    }
    mpfr_clears(s, z, S, zk, term, kln, ks, eps, at, as, (mpfr_ptr)0);
    return out;
}

/* Direct complex power series Li_s(z) = Sum_{k>=1} z^k k^-s, |z| <= 1/2.
 * Returns false if it fails to converge / goes non-finite. */
static bool pl_series_cx(pcx* out, const pcx* s, const pcx* Z, mpfr_prec_t wp) {
    pcx S, zk, ks, term, expo;
    pcx_init(&S, wp); pcx_init(&zk, wp); pcx_init(&ks, wp);
    pcx_init(&term, wp); pcx_init(&expo, wp);
    mpfr_t lnk, eps, at, as;
    mpfr_inits2(wp, lnk, eps, at, as, (mpfr_ptr)0);
    mpfr_set_ui(eps, 1, PRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 8 ? wp - 6 : 1), PRND);

    mpfr_set_ui(S.re, 0, PRND); mpfr_set_ui(S.im, 0, PRND);
    pcx_set(&zk, Z);                              /* z^1 */
    bool ok = true, done = false;

    for (unsigned long k = 1; k < 10000000UL; k++) {
        if (k == 1) {
            pcx_set(&term, &zk);                  /* 1^-s = 1 */
        } else {
            mpfr_set_ui(lnk, k, PRND);
            mpfr_log(lnk, lnk, PRND);
            pcx_mul_r(&expo, s, lnk);              /* s ln k */
            mpfr_neg(expo.re, expo.re, PRND);
            mpfr_neg(expo.im, expo.im, PRND);      /* -s ln k */
            pcx_exp(&ks, &expo, wp);               /* k^-s */
            pcx_mul(&term, &zk, &ks, wp);
        }
        pcx_add(&S, &S, &term);
        if (k >= 3) {
            pcx_abs(at, &term);
            pcx_abs(as, &S);
            mpfr_mul(as, as, eps, PRND);
            if (mpfr_cmp(at, as) < 0) { done = true; break; }
        }
        pcx_mul(&zk, &zk, Z, wp);
    }
    if (done && mpfr_number_p(S.re) && mpfr_number_p(S.im)) pcx_set(out, &S);
    else ok = false;

    pcx_clear(&S); pcx_clear(&zk); pcx_clear(&ks); pcx_clear(&term); pcx_clear(&expo);
    mpfr_clears(lnk, eps, at, as, (mpfr_ptr)0);
    return ok;
}

/* Jonquiere / zeta expansion (DLMF 25.12.11 for non-integer order, 25.12.12 for
 * a positive integer order `intord` >= 2). Valid for |ln z| < 2 Pi. `s` is the
 * (complex) order, used only when intord == 0. Returns false on failure. */
static bool pl_zeta_expansion(pcx* out, const pcx* s, const pcx* Z,
                              long intord, mpfr_prec_t wp) {
    pcx mu, acc, muk, term, w, zk, tmp;
    pcx_init(&mu, wp); pcx_init(&acc, wp); pcx_init(&muk, wp);
    pcx_init(&term, wp); pcx_init(&w, wp); pcx_init(&zk, wp); pcx_init(&tmp, wp);
    mpfr_t fact, eps, at, as;
    mpfr_inits2(wp, fact, eps, at, as, (mpfr_ptr)0);
    mpfr_set_ui(eps, 1, PRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 12 ? wp - 8 : 1), PRND);

    pcx_log(&mu, Z, wp);                          /* mu = ln z */
    bool ok = true;

    /* At z == 1 (mu == 0) the head term degenerates (0 * log(0)); the limit is
     * 0 and the whole expansion collapses to the k = 0 series term zeta(order).
     * Skip the head and let the series produce it. */
    bool mu_zero = mpfr_zero_p(mu.re) && mpfr_zero_p(mu.im);

    if (mu_zero) {
        mpfr_set_zero(acc.re, 1); mpfr_set_zero(acc.im, 1);
    } else if (intord == 0) {
        /* head = Gamma(1 - s) (-mu)^(s-1). */
        pcx oneminus, g, negmu, sm1, p;
        pcx_init(&oneminus, wp); pcx_init(&g, wp); pcx_init(&negmu, wp);
        pcx_init(&sm1, wp); pcx_init(&p, wp);
        mpfr_ui_sub(oneminus.re, 1, s->re, PRND);
        mpfr_neg(oneminus.im, s->im, PRND);
        if (pl_gamma_cx(&g, &oneminus)) {
            mpfr_neg(negmu.re, mu.re, PRND); mpfr_neg(negmu.im, mu.im, PRND);
            mpfr_sub_ui(sm1.re, s->re, 1, PRND); mpfr_set(sm1.im, s->im, PRND);
            pcx_pow(&p, &negmu, &sm1, wp);
            pcx_mul(&acc, &g, &p, wp);
        } else {
            ok = false;
        }
        pcx_clear(&oneminus); pcx_clear(&g); pcx_clear(&negmu);
        pcx_clear(&sm1); pcx_clear(&p);
    } else {
        /* head = mu^(n-1)/(n-1)! (H_{n-1} - ln(-mu)). */
        long nn = intord;
        pcx negmu, lognegmu, coef, munm1;
        pcx_init(&negmu, wp); pcx_init(&lognegmu, wp);
        pcx_init(&coef, wp); pcx_init(&munm1, wp);
        mpfr_t H, t;
        mpfr_inits2(wp, H, t, (mpfr_ptr)0);
        mpfr_set_ui(H, 0, PRND);
        for (long j = 1; j <= nn - 1; j++) {
            mpfr_set_ui(t, 1, PRND);
            mpfr_div_ui(t, t, (unsigned long)j, PRND);
            mpfr_add(H, H, t, PRND);
        }
        mpfr_neg(negmu.re, mu.re, PRND); mpfr_neg(negmu.im, mu.im, PRND);
        pcx_log(&lognegmu, &negmu, wp);
        mpfr_sub(coef.re, H, lognegmu.re, PRND);
        mpfr_neg(coef.im, lognegmu.im, PRND);     /* H - ln(-mu) */
        mpfr_set_ui(munm1.re, 1, PRND); mpfr_set_ui(munm1.im, 0, PRND);
        for (long i = 0; i < nn - 1; i++) { pcx_mul(&tmp, &munm1, &mu, wp); pcx_set(&munm1, &tmp); }
        mpfr_set_ui(t, 1, PRND);                  /* (n-1)! */
        for (long i = 2; i <= nn - 1; i++) mpfr_mul_ui(t, t, (unsigned long)i, PRND);
        pcx_mul(&acc, &munm1, &coef, wp);
        pcx_div_r(&acc, &acc, t);
        mpfr_clears(H, t, (mpfr_ptr)0);
        pcx_clear(&negmu); pcx_clear(&lognegmu); pcx_clear(&coef); pcx_clear(&munm1);
    }

    if (ok) {
        /* + Sum_{k>=0, k != intord-1} zeta(s-k) mu^k/k!. */
        mpfr_set_ui(muk.re, 1, PRND); mpfr_set_ui(muk.im, 0, PRND); /* mu^0 */
        mpfr_set_ui(fact, 1, PRND);                                 /* 0! */
        unsigned long kmax = (unsigned long)(wp * 4 + 200);
        int patience = 0;
        for (unsigned long k = 0; k <= kmax; k++) {
            if (!(intord != 0 && (long)k == intord - 1)) {
                /* w = (order) - k. */
                if (intord != 0) {
                    mpfr_set_si(w.re, intord - (long)k, PRND);
                    mpfr_set_ui(w.im, 0, PRND);
                } else {
                    mpfr_sub_ui(w.re, s->re, k, PRND);
                    mpfr_set(w.im, s->im, PRND);
                }
                if (!pl_zeta_cx(&zk, &w, wp)) { ok = false; break; }
                pcx_mul(&term, &zk, &muk, wp);
                pcx_div_r(&term, &term, fact);
                pcx_add(&acc, &acc, &term);
                pcx_abs(at, &term);
                pcx_abs(as, &acc);
                mpfr_mul(as, as, eps, PRND);
                if (mpfr_cmp(at, as) < 0) {
                    if (++patience >= 2) break;   /* tolerate the zeta(-2j)=0 gaps */
                } else {
                    patience = 0;
                }
            }
            pcx_mul(&tmp, &muk, &mu, wp); pcx_set(&muk, &tmp);
            mpfr_mul_ui(fact, fact, k + 1, PRND);
        }
    }

    if (ok && mpfr_number_p(acc.re) && mpfr_number_p(acc.im)) pcx_set(out, &acc);
    else ok = false;

    mpfr_clears(fact, eps, at, as, (mpfr_ptr)0);
    pcx_clear(&mu); pcx_clear(&acc); pcx_clear(&muk);
    pcx_clear(&term); pcx_clear(&w); pcx_clear(&zk); pcx_clear(&tmp);
    return ok;
}

/* Numeric PolyLog[n, z]. `intord` is the positive integer order (>= 2) or 0.
 * Returns a Real/Complex/MPFR result, or NULL to stay symbolic. */
static Expr* polylog_numeric(Expr* n, Expr* z, long intord) {
    Expr* zero = expr_new_integer(0);
    Expr *nre, *nim, *zre, *zim;
    pl_decompose(n, zero, &nre, &nim);
    pl_decompose(z, zero, &zre, &zim);

    /* All parts must be coercible numbers. */
    double td;
    if (!pl_to_double(nre, &td) || !pl_to_double(zre, &td) ||
        (nim != zero && !pl_to_double(nim, &td)) ||
        (zim != zero && !pl_to_double(zim, &td))) {
        expr_free(zero); return NULL;
    }

    mpfr_prec_t out_prec = pl_out_prec(n, z);
    mpfr_prec_t wp = out_prec + 64;

    bool n_real = (nim == zero);
    bool z_real = (zim == zero);
    double zv = 0.0;
    bool have_zv = z_real && pl_to_double(zre, &zv);

    /* Real fast path: real order, real -1 < z < 1. */
    if (n_real && z_real && have_zv && zv > -1.0 && zv < 1.0) {
        Expr* out = pl_real_series(nre, zre, out_prec, wp);
        expr_free(zero);
        return out;
    }

    /* Branch cut: PolyLog runs from 1 to +Infinity along the real axis. For a
     * real z >= 1 the value is taken continuous from below the cut (Im z -> 0-),
     * matching Mathematica: e.g. Im Li_2(x) = -Pi Log[x] for x > 1. We encode
     * that approach direction as a negative zero imaginary part. */
    bool z_on_cut = z_real && have_zv && zv >= 1.0;

    /* Complex MPFR kernel. */
    pcx s, Z, res;
    pcx_init(&s, wp); pcx_init(&Z, wp); pcx_init(&res, wp);
    bool set = pl_set_mpfr(s.re, nre) && pl_set_mpfr(Z.re, zre);
    if (set) {
        if (nim == zero) mpfr_set_zero(s.im, 1); else set = pl_set_mpfr(s.im, nim);
    }
    if (set) {
        if (zim == zero) mpfr_set_zero(Z.im, z_on_cut ? -1 : 1);
        else set = pl_set_mpfr(Z.im, zim);
    }
    Expr* out = NULL;
    if (set) {
        /* real_only: a real order with a real z < 1 has a real Li_s(z). */
        bool real_only = n_real && z_real && have_zv && zv < 1.0;
        mpfr_t mag, lnmag;
        mpfr_inits2(wp, mag, lnmag, (mpfr_ptr)0);
        pcx_abs(mag, &Z);
        bool got = false;
        if (mpfr_cmp_d(mag, 0.5) <= 0) {
            got = pl_series_cx(&res, &s, &Z, wp);
        } else {
            /* |ln z| < 2 Pi (with a small margin) for the zeta expansion. */
            mpfr_t lz; mpfr_init2(lz, wp);
            {
                pcx muc; pcx_init(&muc, wp);
                pcx_log(&muc, &Z, wp);
                pcx_abs(lz, &muc);
                pcx_clear(&muc);
            }
            if (mpfr_cmp_d(lz, 2.0 * M_PI - 0.05) < 0)
                got = pl_zeta_expansion(&res, &s, &Z, intord, wp);
            mpfr_clear(lz);
        }
        if (got) out = pl_result(res.re, res.im, out_prec, real_only);
        mpfr_clears(mag, lnmag, (mpfr_ptr)0);
    }
    pcx_clear(&s); pcx_clear(&Z); pcx_clear(&res);
    expr_free(zero);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* PolyLog[n, z]                                                       */
/* ------------------------------------------------------------------ */

static Expr* polylog_two_arg(Expr* n, Expr* z) {
    /* PolyLog[n, 0] = 0 (exact zero argument). */
    if (z->type == EXPR_INTEGER && z->data.integer == 0)
        return expr_new_integer(0);

    long m;
    bool isint = pl_integer_order(n, &m);

    /* Order <= 1: rational / logarithmic closed forms (any z, exact or
     * inexact -- the result evaluates numerically when z is inexact). */
    if (isint && m == 1) return pl_order1(z);
    if (isint && m == 0) return pl_order0(z);
    if (isint && m < 0)  {
        Expr* r = pl_neg_order(-m, z);
        if (r) return r;
    }

    /* Order >= 2 with an exact integer order: special argument closed forms. */
    if (isint && m >= 2 && pl_exact_int(n, &m)) {
        long zv;
        if (pl_exact_int(z, &zv) && zv == 1)  return pl_at_one(m);
        if (pl_exact_int(z, &zv) && zv == -1) return pl_at_minus_one(m);
        if (pl_is_half(z)) {
            if (m == 2) return pl_two_half();
            if (m == 3) return pl_three_half();
        }
    }

#ifdef USE_MPFR
    /* Numeric: at least one inexact operand. */
    if (pl_inexact_anywhere(n) || pl_inexact_anywhere(z)) {
        long use = (isint && m >= 2) ? m : 0;
        Expr* out = polylog_numeric(n, z, use);
        if (out) return out;
    }
#else
    (void)isint;
#endif

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

static Expr* polylog_emit_argt(size_t argc) {
    fprintf(stderr,
            "PolyLog::argt: PolyLog called with %zu argument%s; "
            "2 or 3 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_polylog(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 2) return polylog_two_arg(args[0], args[1]);
    if (argc == 3) return NULL;             /* Nielsen S_{n,p}(z): stay symbolic */
    return polylog_emit_argt(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void polylog_init(void) {
    symtab_add_builtin("PolyLog", builtin_polylog);
    symtab_get_def("PolyLog")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

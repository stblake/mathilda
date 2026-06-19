/* Mathilda -- the Lerch transcendent LerchPhi.
 *
 *   LerchPhi[z, s, a]   Phi(z, s, a) = Sum_{k>=0} z^k / (k + a)^s
 *                       (|z| < 1; analytic continuation elsewhere, branch cut
 *                        z in [1, Infinity))
 *
 * For Re(a) < 0 the principal value uses the symmetric power ((k+a)^2)^(-s/2),
 * and any term with k + a = 0 is excluded.  LerchPhi is the common
 * generalization of Zeta, HurwitzZeta and PolyLog:
 *
 *   Phi(1, s, a) = Zeta[s, a],      z Phi(z, s, 1) = PolyLog[s, z].
 *
 * The evaluator routes each kind of argument to the cheapest exact or fastest
 * numeric path (mirroring src/special_functions/{zeta,hurwitzzeta,polylog}.c):
 *
 *   exact reductions (any z, s, a):
 *       z = 0                  ->  a^-s                        (the k = 0 term)
 *       s = 0                  ->  1/(1 - z)                   (geometric sum)
 *       z = 1                  ->  Zeta[s, a]
 *       z = -1                 ->  2^-s (Zeta[s,a/2] - Zeta[s,(a+1)/2])
 *       a positive integer m   ->  z^-m (PolyLog[s,z] - Sum_{j<m} z^j j^-s)
 *       s negative integer -n  ->  (z d/dz + a)^n [1/(1-z)]    (rational in z)
 *   options:
 *       IncludeSingularTerm->True at a non-positive integer a  ->  ComplexInfinity
 *       DoublyInfinite->True   ->  Phi(z,s,a) + z^-1 Phi(1/z, s, 1-a)
 *   numeric (>= 1 inexact operand):
 *       |z| < 1                ->  complex-MPFR power series
 *       |z| > 1                ->  diverges, stays symbolic (no continuation)
 *   everything else            ->  stays symbolic (return NULL)
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "lerchphi.h"
#include "sym_names.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "arithmetic.h"   /* is_rational, make_rational, is_complex, make_complex */
#include "numeric.h"      /* numeric_min_inexact_bits */
#include "print.h"        /* expr_to_string (nonopt diagnostic) */
#include "attr.h"
#include "eval.h"          /* eval_and_free */
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#define LRND MPFR_RNDN
#endif

/* M_PI is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Largest positive integer a expanded into the PolyLog reduction. */
#define LP_POS_A_CAP   100000L
/* Largest n for which LerchPhi[z, -n, a] is expanded into its rational form. */
#define LP_NEG_S_CAP   200L

/* ------------------------------------------------------------------ */
/* Small predicates / coercions                                        */
/* ------------------------------------------------------------------ */

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool lp_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if any leaf of `e` (including the parts of a Complex) is inexact. */
static bool lp_has_inexact(Expr* e) {
    if (!e) return false;
    if (lp_is_inexact(e)) return true;
    Expr *re, *im;
    if (is_complex(e, &re, &im)) return lp_is_inexact(re) || lp_is_inexact(im);
    return false;
}

/* Extract an exact machine integer from `e` (Integer, or BigInt fitting a
 * signed long). Inexact reals do NOT count. */
static bool lp_exact_int(const Expr* e, long* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (long)e->data.integer; return true; }
    if (e->type == EXPR_BIGINT && mpz_fits_slong_p(e->data.bigint)) {
        *out = mpz_get_si(e->data.bigint); return true;
    }
    return false;
}

/* True if `e` is exactly the integer value `v`. */
static bool lp_is_int(const Expr* e, long v) {
    long n;
    return lp_exact_int(e, &n) && n == v;
}

/* ------------------------------------------------------------------ */
/* Small Expr builders                                                 */
/* ------------------------------------------------------------------ */

static Expr* lp_pow(Expr* base, Expr* expo) {     /* base^expo, both consumed */
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ base, expo }, 2);
}
static Expr* lp_times(Expr* a, Expr* b) {         /* a b, both consumed */
    return expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ a, b }, 2);
}
static Expr* lp_neg_s(Expr* s) {                  /* -s = Times[-1, s], s consumed */
    return lp_times(expr_new_integer(-1), s);
}
static Expr* lp_zeta2(Expr* s, Expr* a) {         /* Zeta[s, a], both consumed */
    return expr_new_function(expr_new_symbol(SYM_Zeta), (Expr*[]){ s, a }, 2);
}

/* ------------------------------------------------------------------ */
/* Exact reductions                                                    */
/* ------------------------------------------------------------------ */

/* z = 0: only the k = 0 term survives, Phi(0, s, a) = a^-s. */
static Expr* lp_z_zero(Expr* s, Expr* a) {
    return eval_and_free(lp_pow(expr_copy(a), lp_neg_s(expr_copy(s))));
}

/* s = 0: Phi(z, 0, a) = Sum z^k = 1/(1 - z), independent of a. */
static Expr* lp_s_zero(Expr* z) {
    Expr* omz = expr_new_function(expr_new_symbol(SYM_Subtract),
                    (Expr*[]){ expr_new_integer(1), expr_copy(z) }, 2);
    return eval_and_free(lp_pow(omz, expr_new_integer(-1)));
}

/* z = 1: Phi(1, s, a) = Zeta[s, a] (same symmetric convention as two-arg Zeta). */
static Expr* lp_z_one(Expr* s, Expr* a) {
    return eval_and_free(lp_zeta2(expr_copy(s), expr_copy(a)));
}

/* z = -1: Phi(-1, s, a) = 2^-s (Zeta[s, a/2] - Zeta[s, (a+1)/2]). */
static Expr* lp_z_minus_one(Expr* s, Expr* a) {
    Expr* twonegs = lp_pow(expr_new_integer(2), lp_neg_s(expr_copy(s)));
    Expr* ah  = lp_times(expr_copy(a), make_rational(1, 2));               /* a/2 */
    Expr* ap1 = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(a), expr_new_integer(1) }, 2);
    Expr* ah2 = lp_times(ap1, make_rational(1, 2));                        /* (a+1)/2 */
    Expr* z1  = lp_zeta2(expr_copy(s), ah);
    Expr* z2  = lp_zeta2(expr_copy(s), ah2);
    Expr* diff = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ z1, lp_times(expr_new_integer(-1), z2) }, 2);
    return eval_and_free(lp_times(twonegs, diff));
}

/* a = m, a positive integer: shift the index down to the PolyLog series,
 *   Phi(z, s, m) = z^-m (PolyLog[s, z] - Sum_{j=1}^{m-1} z^j j^-s).
 * Reuses the PolyLog engine (so exact / symbolic / numeric all follow). NULL
 * past the cap. */
static Expr* lp_pos_int_a(Expr* z, Expr* s, long m) {
    if (m < 1 || m > LP_POS_A_CAP) return NULL;

    Expr* poly = expr_new_function(expr_new_symbol(SYM_PolyLog),
                     (Expr*[]){ expr_copy(s), expr_copy(z) }, 2);
    Expr* bracket;
    if (m == 1) {
        bracket = poly;                            /* empty subtracted sum */
    } else {
        size_t terms = (size_t)m;                  /* poly + (m-1) subtracted terms */
        Expr** parts = (Expr**)malloc(terms * sizeof(Expr*));
        if (!parts) { expr_free(poly); return NULL; }
        parts[0] = poly;
        for (long j = 1; j < m; j++) {
            /* -(z^j j^-s) = Times[-1, Power[z, j], Power[j, -s]]. */
            Expr* zj = lp_pow(expr_copy(z), expr_new_integer(j));
            Expr* js = lp_pow(expr_new_integer(j), lp_neg_s(expr_copy(s)));
            parts[(size_t)j] = expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ expr_new_integer(-1), zj, js }, 3);
        }
        bracket = expr_new_function(expr_new_symbol(SYM_Plus), parts, terms);
        free(parts);
    }
    Expr* zinv = lp_pow(expr_copy(z), expr_new_integer(-m));
    return eval_and_free(lp_times(zinv, bracket));
}

/* s = -n, n >= 1: Phi(z, -n, a) = Sum_{k>=0} (k+a)^n z^k = (z d/dz + a)^n
 * applied to 1/(1-z), a rational function of z.  Built by applying the Euler
 * operator theta = z d/dz n times, then Together. NULL past the cap. */
static Expr* lp_neg_int_s(Expr* z, Expr* a, long n) {
    if (n < 1 || n > LP_NEG_S_CAP) return NULL;

    /* Apply the operator against a fresh placeholder symbol so that
     * differentiation works even when z is a concrete number (a literal value
     * would collapse 1/(1-z) to a constant and lose every derivative); the real
     * z is substituted only after the operator has been applied. */
    Expr* zv = expr_new_symbol("$LerchPhiZ");

    /* cur = 1/(1 - zv). */
    Expr* omz = expr_new_function(expr_new_symbol(SYM_Subtract),
                    (Expr*[]){ expr_new_integer(1), expr_copy(zv) }, 2);
    Expr* cur = eval_and_free(lp_pow(omz, expr_new_integer(-1)));

    for (long i = 0; i < n; i++) {
        /* dz = D[cur, zv]. */
        Expr* dz = eval_and_free(expr_new_function(expr_new_symbol(SYM_D),
                       (Expr*[]){ expr_copy(cur), expr_copy(zv) }, 2));
        Expr* term1 = lp_times(expr_copy(zv), dz);    /* zv d/dzv cur */
        Expr* term2 = lp_times(expr_copy(a), cur);    /* a cur (consumes cur) */
        cur = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                       (Expr*[]){ term1, term2 }, 2));
    }

    /* Substitute zv -> z, then Together. */
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                     (Expr*[]){ expr_copy(zv), expr_copy(z) }, 2);
    cur = eval_and_free(expr_new_function(expr_new_symbol("ReplaceAll"),
                     (Expr*[]){ cur, rule }, 2));
    expr_free(zv);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Together),
                       (Expr*[]){ cur }, 1));
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ */
/* Complex-MPFR toolkit (pairs of mpfr_t; no MPC library available).   */
/* Mirrors the `pcx`/`zcx` helpers -- alias-safe, explicit precision p. */
/* ------------------------------------------------------------------ */

typedef struct { mpfr_t re, im; } lcx;

static void lcx_init(lcx* z, mpfr_prec_t p) { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void lcx_clear(lcx* z)               { mpfr_clear(z->re);    mpfr_clear(z->im);    }
static void lcx_set(lcx* d, const lcx* s)   { mpfr_set(d->re, s->re, LRND); mpfr_set(d->im, s->im, LRND); }

static void lcx_add(lcx* out, const lcx* a, const lcx* b) {
    mpfr_add(out->re, a->re, b->re, LRND);
    mpfr_add(out->im, a->im, b->im, LRND);
}

static void lcx_abs(mpfr_t mag, const lcx* z) { mpfr_hypot(mag, z->re, z->im, LRND); }

static void lcx_mul(lcx* out, const lcx* a, const lcx* b, mpfr_prec_t p) {
    mpfr_t ac, bd, ad, bc;
    mpfr_inits2(p, ac, bd, ad, bc, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, LRND);
    mpfr_mul(bd, a->im, b->im, LRND);
    mpfr_mul(ad, a->re, b->im, LRND);
    mpfr_mul(bc, a->im, b->re, LRND);
    mpfr_sub(out->re, ac, bd, LRND);
    mpfr_add(out->im, ad, bc, LRND);
    mpfr_clears(ac, bd, ad, bc, (mpfr_ptr)0);
}

static void lcx_exp(lcx* out, const lcx* z, mpfr_prec_t p) {
    mpfr_t ea, c, s;
    mpfr_inits2(p, ea, c, s, (mpfr_ptr)0);
    mpfr_exp(ea, z->re, LRND);
    mpfr_sin_cos(s, c, z->im, LRND);
    mpfr_mul(out->re, ea, c, LRND);
    mpfr_mul(out->im, ea, s, LRND);
    mpfr_clears(ea, c, s, (mpfr_ptr)0);
}

/* out = Log(z), principal branch: ln|z| + i Arg(z). */
static void lcx_log(lcx* out, const lcx* z, mpfr_prec_t p) {
    mpfr_t mag, ang;
    mpfr_inits2(p, mag, ang, (mpfr_ptr)0);
    mpfr_hypot(mag, z->re, z->im, LRND);
    mpfr_atan2(ang, z->im, z->re, LRND);
    mpfr_log(out->re, mag, LRND);
    mpfr_set(out->im, ang, LRND);
    mpfr_clears(mag, ang, (mpfr_ptr)0);
}

/* out = base^expo (principal branch) = exp(expo * Log(base)). */
static void lcx_pow(lcx* out, const lcx* base, const lcx* expo, mpfr_prec_t p) {
    lcx lg, prod;
    lcx_init(&lg, p); lcx_init(&prod, p);
    lcx_log(&lg, base, p);
    lcx_mul(&prod, expo, &lg, p);
    lcx_exp(out, &prod, p);
    lcx_clear(&lg); lcx_clear(&prod);
}

/* ------------------------------------------------------------------ */
/* MPFR coercion / result builders                                     */
/* ------------------------------------------------------------------ */

/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool lp_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, LRND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          LRND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        LRND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          LRND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, LRND);
        mpfr_div_si(out, out, (long)d, LRND);
        return true;
    }
    return false;
}

/* Result precision (bits) under Mathematica contagion: the minimum precision
 * among the inexact leaves of z, s, a, floored at machine (53). */
static mpfr_prec_t lp_out_prec(const Expr* z, const Expr* s, const Expr* a) {
    long p[3];
    p[0] = numeric_min_inexact_bits((Expr*)z);
    p[1] = numeric_min_inexact_bits((Expr*)s);
    p[2] = numeric_min_inexact_bits((Expr*)a);
    long m = 0;
    for (int i = 0; i < 3; i++)
        if (p[i] && (m == 0 || p[i] < m)) m = p[i];
    return m < 53 ? 53 : (mpfr_prec_t)m;
}

/* Build the numeric result: Real / Complex[Real,Real] at machine precision,
 * MPFR parts at arbitrary precision; `real_only` forces a real result. */
static Expr* lp_result(const mpfr_t re, const mpfr_t im,
                       mpfr_prec_t out_prec, bool real_only) {
    if (out_prec <= 53) {
        Expr* rr = expr_new_real(mpfr_get_d(re, LRND));
        if (real_only) return rr;
        return make_complex(rr, expr_new_real(mpfr_get_d(im, LRND)));
    }
    Expr* rr = expr_new_mpfr_bits(out_prec);
    mpfr_set(rr->data.mpfr, re, LRND);
    if (real_only) return rr;
    Expr* ii = expr_new_mpfr_bits(out_prec);
    mpfr_set(ii->data.mpfr, im, LRND);
    return make_complex(rr, ii);
}

/* Decompose a numeric value into borrowed (re, im); plain real -> (value, zero). */
static void lp_decompose(Expr* v, Expr* zero, Expr** re, Expr** im) {
    Expr *r, *i;
    if (is_complex(v, &r, &i)) { *re = r; *im = i; }
    else { *re = v; *im = zero; }
}

/* ------------------------------------------------------------------ */
/* Numeric power-series kernel                                         */
/* ------------------------------------------------------------------ */

enum { LP_OK = 0, LP_NOCONV = 1, LP_SINGULAR = 2 };

/* Phi(z, s, a) = Sum_{k>=0} z^k (k+a)^-s into `out` (precision wp).  Uses the
 * symmetric power ((k+a)^2)^(-s/2) when Re(a) < 0, else the principal (k+a)^-s.
 * A term with k + a == 0 is excluded unless `include_singular` (then the result
 * is singular).  Converges only for |z| < 1 (the caller guards |z| > 1).
 * Returns LP_OK / LP_NOCONV / LP_SINGULAR. */
static int lp_series_cx(lcx* out, const lcx* Z, const lcx* s, const lcx* a,
                        bool include_singular, mpfr_prec_t wp) {
    bool symmetric = (mpfr_sgn(a->re) < 0);

    lcx S, zk, base, b2, val, term, expo;
    lcx_init(&S, wp); lcx_init(&zk, wp); lcx_init(&base, wp); lcx_init(&b2, wp);
    lcx_init(&val, wp); lcx_init(&term, wp); lcx_init(&expo, wp);

    /* expo = -s (principal) or -s/2 (symmetric). */
    mpfr_neg(expo.re, s->re, LRND);
    mpfr_neg(expo.im, s->im, LRND);
    if (symmetric) { mpfr_div_2ui(expo.re, expo.re, 1, LRND);
                     mpfr_div_2ui(expo.im, expo.im, 1, LRND); }

    mpfr_t eps, mag, at, as;
    mpfr_inits2(wp, eps, mag, at, as, (mpfr_ptr)0);
    mpfr_set_ui(eps, 1, LRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 12 ? wp - 8 : 1), LRND);

    mpfr_set_ui(S.re, 0, LRND); mpfr_set_ui(S.im, 0, LRND);
    mpfr_set_ui(zk.re, 1, LRND); mpfr_set_ui(zk.im, 0, LRND);   /* z^0 */

    int code = LP_NOCONV;
    int patience = 0;
    unsigned long hard_cap = (unsigned long)wp * 64UL + 100000UL;
    for (unsigned long k = 0; k < hard_cap; k++) {
        mpfr_add_ui(base.re, a->re, k, LRND);
        mpfr_set(base.im, a->im, LRND);
        lcx_abs(mag, &base);
        if (mpfr_sgn(mag) == 0) {                  /* k + a == 0 */
            if (include_singular) { code = LP_SINGULAR; break; }
            lcx_mul(&zk, &zk, Z, wp);              /* skip term, advance z^k */
            continue;
        }
        if (symmetric) {
            lcx_mul(&b2, &base, &base, wp);        /* (k+a)^2 */
            lcx_pow(&val, &b2, &expo, wp);         /* ((k+a)^2)^(-s/2) */
        } else {
            lcx_pow(&val, &base, &expo, wp);       /* (k+a)^-s */
        }
        lcx_mul(&term, &zk, &val, wp);
        lcx_add(&S, &S, &term);

        if (k >= 2) {
            lcx_abs(at, &term);
            lcx_abs(as, &S);
            mpfr_mul(as, as, eps, LRND);
            if (mpfr_cmp(at, as) < 0) {
                if (++patience >= 2) { code = LP_OK; break; }
            } else {
                patience = 0;
            }
        }
        lcx_mul(&zk, &zk, Z, wp);
    }

    if (code == LP_OK && !(mpfr_number_p(S.re) && mpfr_number_p(S.im)))
        code = LP_NOCONV;
    if (code == LP_OK) lcx_set(out, &S);

    mpfr_clears(eps, mag, at, as, (mpfr_ptr)0);
    lcx_clear(&S); lcx_clear(&zk); lcx_clear(&base); lcx_clear(&b2);
    lcx_clear(&val); lcx_clear(&term); lcx_clear(&expo);
    return code;
}

/* Numeric LerchPhi[z, s, a] (at least one inexact operand).  Returns a
 * Real/Complex/MPFR result, ComplexInfinity at an included singular term, or
 * NULL to stay symbolic (non-numeric parts, or |z| > 1 where the series
 * diverges and no continuation is implemented). */
static Expr* lp_numeric(Expr* z, Expr* s, Expr* a, bool include_singular) {
    Expr* zero = expr_new_integer(0);
    Expr *zre, *zim, *sre, *sim, *are, *aim;
    lp_decompose(z, zero, &zre, &zim);
    lp_decompose(s, zero, &sre, &sim);
    lp_decompose(a, zero, &are, &aim);

    mpfr_prec_t out_prec = lp_out_prec(z, s, a);
    mpfr_prec_t wp = out_prec + 64;

    bool real_only = (zim == zero) && (sim == zero) && (aim == zero);

    lcx Z, sc, ac, res;
    lcx_init(&Z, wp); lcx_init(&sc, wp); lcx_init(&ac, wp); lcx_init(&res, wp);
    bool ok = lp_set_mpfr(Z.re, zre) && lp_set_mpfr(Z.im, zim) &&
              lp_set_mpfr(sc.re, sre) && lp_set_mpfr(sc.im, sim) &&
              lp_set_mpfr(ac.re, are) && lp_set_mpfr(ac.im, aim);

    Expr* out = NULL;
    if (ok) {
        /* |z| > 1: the defining series diverges and there is no analytic
         * continuation here -- leave symbolic.  (|z| = 1 is attempted; it
         * converges when Re(s) > 1.) */
        mpfr_t zmag; mpfr_init2(zmag, wp);
        lcx_abs(zmag, &Z);
        bool diverges = mpfr_cmp_d(zmag, 1.0 + 1e-12) > 0;
        mpfr_clear(zmag);
        if (!diverges) {
            int code = lp_series_cx(&res, &Z, &sc, &ac, include_singular, wp);
            if (code == LP_OK)            out = lp_result(res.re, res.im, out_prec, real_only);
            else if (code == LP_SINGULAR) out = expr_new_symbol(SYM_ComplexInfinity);
        }
    }
    lcx_clear(&Z); lcx_clear(&sc); lcx_clear(&ac); lcx_clear(&res);
    expr_free(zero);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Core dispatch                                                       */
/* ------------------------------------------------------------------ */

static Expr* lerchphi_core(Expr* z, Expr* s, Expr* a,
                           bool doubly_infinite, bool include_singular);

/* DoublyInfinite -> True:  Phi_DI(z,s,a) = Phi(z,s,a) + z^-1 Phi(1/z, s, 1-a).
 * The two single sums re-enter the builtin with default options. */
static Expr* lp_doubly_infinite(Expr* z, Expr* s, Expr* a, bool include_singular) {
    Expr* main = lerchphi_core(z, s, a, false, include_singular);
    if (!main) return NULL;                        /* one half symbolic -> stay symbolic */

    Expr* zinv = lp_pow(expr_copy(z), expr_new_integer(-1));        /* 1/z */
    Expr* oma  = expr_new_function(expr_new_symbol(SYM_Subtract),
                     (Expr*[]){ expr_new_integer(1), expr_copy(a) }, 2);  /* 1 - a */
    Expr* refl = expr_new_function(expr_new_symbol(SYM_LerchPhi),
                     (Expr*[]){ expr_copy(zinv), expr_copy(s), oma }, 3);
    Expr* tail = eval_and_free(lp_times(zinv, refl));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ main, tail }, 2));
}

/* Single-sum dispatch (DoublyInfinite already peeled off). */
static Expr* lerchphi_single(Expr* z, Expr* s, Expr* a, bool include_singular) {
    long m;

    /* IncludeSingularTerm at a non-positive integer a -> the k = -a term blows
     * up (0^-s) for any s. */
    if (include_singular && lp_exact_int(a, &m) && m <= 0)
        return expr_new_symbol(SYM_ComplexInfinity);

    /* z = 0: a^-s. */
    if (lp_is_int(z, 0)) return lp_z_zero(s, a);

    /* s = 0: 1/(1 - z) (also handles z = 1 -> ComplexInfinity). */
    if (lp_is_int(s, 0)) return lp_s_zero(z);

    /* z = 1: Zeta[s, a]. */
    if (lp_is_int(z, 1)) return lp_z_one(s, a);

    /* z = -1: 2^-s (Zeta[s,a/2] - Zeta[s,(a+1)/2]). */
    if (lp_is_int(z, -1)) return lp_z_minus_one(s, a);

    /* a positive integer: PolyLog reduction. */
    if (lp_exact_int(a, &m) && m >= 1) {
        Expr* r = lp_pos_int_a(z, s, m);
        if (r) return r;
    }

    /* s negative integer: rational function of z. */
    if (lp_exact_int(s, &m) && m < 0) {
        Expr* r = lp_neg_int_s(z, a, -m);
        if (r) return r;
    }

#ifdef USE_MPFR
    /* Numeric: at least one inexact operand. */
    if (lp_has_inexact(z) || lp_has_inexact(s) || lp_has_inexact(a)) {
        Expr* out = lp_numeric(z, s, a, include_singular);
        if (out) return out;
    }
#else
    (void)include_singular;
#endif

    return NULL; /* leave symbolic */
}

static Expr* lerchphi_core(Expr* z, Expr* s, Expr* a,
                           bool doubly_infinite, bool include_singular) {
    if (doubly_infinite) return lp_doubly_infinite(z, s, a, include_singular);
    return lerchphi_single(z, s, a, include_singular);
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

static Expr* lp_emit_argrx(size_t argc) {
    fprintf(stderr,
            "LerchPhi::argrx: LerchPhi called with %zu argument%s; "
            "3 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Mirror Wolfram's nonopt surface text for a bad option beyond position 3. */
static Expr* lp_emit_nonopt(Expr* bad, Expr* res) {
    char* bad_str  = expr_to_string(bad);
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "LerchPhi::nonopt: Options expected (instead of %s) beyond position 3 "
            "in %s. An option must be a rule or a list of rules.\n",
            bad_str ? bad_str : "?", call_str ? call_str : "?");
    free(bad_str);
    free(call_str);
    return NULL;
}

/* True if `opt` is Rule[symbol, value]; fills *name and *val when so. */
static bool lp_as_option(Expr* opt, const char** name, Expr** val) {
    if (opt->type == EXPR_FUNCTION &&
        opt->data.function.head->type == EXPR_SYMBOL &&
        opt->data.function.head->data.symbol == SYM_Rule &&
        opt->data.function.arg_count == 2 &&
        opt->data.function.args[0]->type == EXPR_SYMBOL) {
        *name = opt->data.function.args[0]->data.symbol;
        *val  = opt->data.function.args[1];
        return true;
    }
    return false;
}

Expr* builtin_lerchphi(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc < 3) return lp_emit_argrx(argc);

    bool doubly_infinite = false, include_singular = false;
    Expr* last_bad = NULL;
    for (size_t i = 3; i < argc; i++) {
        const char* name; Expr* val;
        if (!lp_as_option(args[i], &name, &val)) { last_bad = args[i]; continue; }
        bool is_true = (val->type == EXPR_SYMBOL && val->data.symbol == SYM_True);
        if (name == SYM_DoublyInfinite)            doubly_infinite = is_true;
        else if (name == SYM_IncludeSingularTerm)  include_singular = is_true;
        else                                       last_bad = args[i];
    }
    if (last_bad) return lp_emit_nonopt(last_bad, res);

    return lerchphi_core(args[0], args[1], args[2], doubly_infinite, include_singular);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void lerchphi_init(void) {
    symtab_add_builtin("LerchPhi", builtin_lerchphi);
    symtab_get_def("LerchPhi")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}

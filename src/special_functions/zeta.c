/* Mathilda -- the Riemann and Hurwitz zeta functions.
 *
 *   Zeta[s]      Riemann zeta      zeta(s)   = Sum_{k>=1} k^-s          (Re s > 1)
 *   Zeta[s, a]   Hurwitz zeta      zeta(s,a) = Sum_{k>=0} (k+a)^-s      (Re s > 1)
 *
 * Both are defined elsewhere by analytic continuation; the evaluator routes
 * each kind of argument to the cheapest exact or fastest numeric path:
 *
 *   exact integer s         ->  closed form:
 *                                 s = 1        : ComplexInfinity (pole)
 *                                 s = 0        : -1/2
 *                                 s = 2n > 0   : rational * Pi^(2n)   (Bernoulli)
 *                                 s = -m < 0   : rational            (Bernoulli)
 *                                 s = 2n+1 > 0 : stays symbolic (no closed form)
 *   exact Hurwitz, integer a -> Zeta[s] - Sum_{k=1}^{a-1} k^-s
 *   exact Hurwitz, a <= 0    -> Sum_{k+a<0} |k+a|^-s + Zeta[s, a+m]
 *                               (the symmetric convention; see below)
 *   Zeta[0, a], rational a>0 -> 1/2 - a
 *   machine / MPFR real s    -> MPFR mpfr_zeta (Riemann only)
 *   complex s, or any a != 1 -> Euler-Maclaurin complex-MPFR kernel
 *   everything else          -> stays symbolic (return NULL)
 *
 * MPFR provides mpfr_zeta for real Riemann zeta only -- it has no Hurwitz and
 * no complex zeta -- so the Hurwitz / complex kernel is implemented here from
 * the Euler-Maclaurin summation formula (DLMF 25.11.5):
 *
 *   zeta(s,a) = Sum_{k=0}^{N-1} (a+k)^-s
 *             + (a+N)^(1-s)/(s-1)
 *             + 1/2 (a+N)^-s
 *             + Sum_{j>=1} B_{2j}/(2j)! (s)_{2j-1} (a+N)^(-s-2j+1)
 *
 * with (s)_{2j-1} the rising factorial. N is chosen from the working precision
 * and |s|; the correction series is truncated at its optimal (smallest) term.
 *
 * For a <= 0 the defining sum is taken in the symmetric form
 * Sum_k ((k+a)^2)^(-s/2), and the k + a = 0 term is dropped unless
 * IncludeSingularTerm -> True, which makes Zeta[s, a] ComplexInfinity at a
 * non-positive integer a whenever Re s > 0. (At s = 0 that term is 0^0 = 1 and
 * is always counted; for Re s < 0 it is 0.)
 *
 * Options: IncludeSingularTerm -> False.
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "zeta.h"
#include "sym_names.h"
#include "interval.h"

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
#include "flint_num_bridge.h" /* flint_num_zeta (arb acb_dirichlet_zeta) */
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "options.h"   /* options_parse_bools */
#include "ndarray.h"   /* is_ndarray, ndarray_delist_and_reeval */
#include "sf_scalar.h"       /* sf_machine_scalar_args: the machine-scalar gate */
#include "sf_machine_zeta.h" /* sf_zeta_r/_c, sf_zeta2_r/_c: the shared kernels */
#include "deriv.h"           /* deriv_register_numeric */

#ifdef USE_MPFR
#include <mpfr.h>
#define ZRND MPFR_RNDN
#endif

#ifdef USE_FLINT
#include <flint/fmpq.h>
#include <flint/arith.h>       /* arith_bernoulli_number */
#endif

/* M_PI is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Largest |s| (exact integer) for which we build the Bernoulli closed form;
 * beyond this the Bernoulli numbers get unwieldy, so stay symbolic / numeric. */
#define ZETA_EXACT_INT_CAP 10000L
/* Largest positive integer a for which exact Hurwitz expands to a partial sum. */
#define ZETA_HURWITZ_A_CAP 100000L

/* ------------------------------------------------------------------ */
/* Small predicates / coercions                                        */
/* ------------------------------------------------------------------ */

/* True if `e` is exactly the symbol `name`. */
static bool zeta_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool zeta_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* Extract an exact machine integer from `e` (Integer, or BigInt fitting a
 * signed long). Inexact reals do NOT count -- the exact path is only for
 * genuinely exact integer arguments. */
static bool zeta_exact_int(const Expr* e, long* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (long)e->data.integer; return true; }
    if (e->type == EXPR_BIGINT && mpz_fits_slong_p(e->data.bigint)) {
        *out = mpz_get_si(e->data.bigint); return true;
    }
    return false;
}

/* True if `e` is exactly the rational 1/2. */
static bool zeta_is_half(const Expr* e) {
    int64_t n, d;
    return is_rational(e, &n, &d) && n == 1 && d == 2;
}

/* ------------------------------------------------------------------ */
/* Bernoulli numbers B_k (exact, lazily cached, process-lifetime)      */
/*                                                                     */
/* Recurrence: B_0 = 1, and for m >= 1                                 */
/*   B_m = -1/(m+1) Sum_{k=0}^{m-1} C(m+1, k) B_k.                      */
/* (Same construction as src/polygamma.c; replicated to stay            */
/* self-contained -- the cache there is file-static.)                  */
/* ------------------------------------------------------------------ */

static mpq_t* g_zbern = NULL;
static size_t g_zbern_len = 0;

static void zbern_ensure(size_t upto) {
    if (upto + 1 <= g_zbern_len) return;
    size_t newlen = upto + 1;
    mpq_t* grown = (mpq_t*)realloc(g_zbern, newlen * sizeof(mpq_t));
    if (!grown) return; /* out of memory: leave cache as-is, callers degrade */
    g_zbern = grown;

    mpq_t sum, term, factor;
    mpq_inits(sum, term, factor, (mpq_ptr)0);
    mpz_t binz;
    mpz_init(binz);

    for (size_t m = g_zbern_len; m < newlen; m++) {
        mpq_init(g_zbern[m]);
        if (m == 0) { mpq_set_ui(g_zbern[m], 1, 1); continue; }
        mpq_set_ui(sum, 0, 1);
        for (size_t k = 0; k < m; k++) {
            mpz_bin_uiui(binz, (unsigned long)(m + 1), (unsigned long)k);
            mpq_set_z(factor, binz);
            mpq_mul(term, factor, g_zbern[k]);
            mpq_add(sum, sum, term);
        }
        mpq_set_si(term, -1, (unsigned long)(m + 1));
        mpq_canonicalize(term);
        mpq_mul(g_zbern[m], sum, term);
        mpq_canonicalize(g_zbern[m]);
    }
    g_zbern_len = newlen;
    mpz_clear(binz);
    mpq_clears(sum, term, factor, (mpq_ptr)0);
}

#ifdef USE_FLINT
/* With FLINT each requested B_idx is computed directly (arith_bernoulli_number:
 * tabulated, zeta-based or multi-modular by size) instead of by the O(n^2)
 * recurrence over every smaller index, and kept in a sparse cache: Zeta[400]
 * cold went from 31 ms to well under 1 ms. */
static mpq_t* g_zbq = NULL;
static unsigned char* g_zbq_set = NULL;
static size_t g_zbq_len = 0;

static bool zbern_flint(mpq_t out, size_t idx) {
    if (idx >= g_zbq_len) {
        size_t newlen = idx + 1 > 2 * g_zbq_len ? idx + 1 : 2 * g_zbq_len;
        mpq_t* q = (mpq_t*)realloc(g_zbq, newlen * sizeof(mpq_t));
        if (!q) return false;
        g_zbq = q;
        unsigned char* f = (unsigned char*)realloc(g_zbq_set, newlen);
        if (!f) return false;
        g_zbq_set = f;
        for (size_t i = g_zbq_len; i < newlen; i++) { mpq_init(g_zbq[i]); g_zbq_set[i] = 0; }
        g_zbq_len = newlen;
    }
    if (!g_zbq_set[idx]) {
        fmpq_t b;
        fmpq_init(b);
        arith_bernoulli_number(b, (ulong)idx);
        fmpq_get_mpq(g_zbq[idx], b);
        fmpq_clear(b);
        g_zbq_set[idx] = 1;
    }
    mpq_set(out, g_zbq[idx]);
    return true;
}
#endif

/* B_idx into `out`. Odd indices above 1 are exactly zero. */
static void zbern_get_q(mpq_t out, size_t idx) {
    if (idx > 1 && (idx & 1u)) { mpq_set_ui(out, 0, 1); return; }
#ifdef USE_FLINT
    if (zbern_flint(out, idx)) return;
#endif
    zbern_ensure(idx);
    if (idx < g_zbern_len) mpq_set(out, g_zbern[idx]);
    else mpq_set_ui(out, 0, 1);
}

/* ------------------------------------------------------------------ */
/* Exact-value builders                                                */
/* ------------------------------------------------------------------ */

/* n! into out. */
static void zeta_factorial(mpz_t out, long n) {
    mpz_fac_ui(out, (unsigned long)n);
}

/* A canonical Integer/BigInt/Rational from a canonical mpq. */
static Expr* zeta_expr_from_mpq(const mpq_t q) {
    mpz_t num, den;
    mpz_init(num); mpz_init(den);
    mpq_get_num(num, q);
    mpq_get_den(den, q);
    Expr* out = make_rational_mpz(num, den);   /* q is canonical: no gcd work */
    mpz_clear(num); mpz_clear(den);
    return out;
}

/* zeta(2n) = (-1)^(n+1) B_{2n} 2^(2n-1) Pi^(2n) / (2n)!,  for n >= 1.
 * Returns rational * Pi^(2n), or NULL if 2n exceeds the cap. */
/* The built zeta(2n) values, by n: an expression is immutable and expr_copy
 * is a reference bump, so a repeated Zeta[2n] costs no rebuilding of
 * rational * Pi^(2n). Process-lifetime, like the Bernoulli cache. */
static Expr** g_zeven = NULL;
static size_t g_zeven_len = 0;

static Expr* zeta_exact_even_build(long two_n);

static Expr* zeta_exact_even(long two_n) {
    if (two_n > ZETA_EXACT_INT_CAP) return NULL;
    size_t n = (size_t)(two_n / 2);
    if (n >= g_zeven_len) {
        size_t newlen = n + 1 > 2 * g_zeven_len ? n + 1 : 2 * g_zeven_len;
        Expr** grown = (Expr**)realloc(g_zeven, newlen * sizeof(Expr*));
        if (!grown) return zeta_exact_even_build(two_n);
        for (size_t i = g_zeven_len; i < newlen; i++) grown[i] = NULL;
        g_zeven = grown;
        g_zeven_len = newlen;
    }
    if (!g_zeven[n]) g_zeven[n] = zeta_exact_even_build(two_n);
    return g_zeven[n] ? expr_copy(g_zeven[n]) : NULL;
}

static Expr* zeta_exact_even_build(long two_n) {
    long n = two_n / 2;

    mpq_t c, tmp;
    mpz_t f2n;
    mpq_inits(c, tmp, (mpq_ptr)0);
    mpz_init(f2n);

    zbern_get_q(c, (size_t)two_n);                 /* B_{2n} */
    mpq_mul_2exp(c, c, (mp_bitcnt_t)(two_n - 1));  /* 2^(2n-1) */
    zeta_factorial(f2n, two_n);                    /* (2n)! */
    mpq_set_z(tmp, f2n);
    mpq_div(c, c, tmp);
    if ((n % 2) == 0) mpq_neg(c, c);               /* (-1)^(n+1) */
    mpq_canonicalize(c);

    Expr* crat = zeta_expr_from_mpq(c);
    Expr* pipow = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_new_symbol(SYM_Pi), expr_new_integer(two_n) }, 2);
    Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ crat, pipow }, 2));

    mpq_clears(c, tmp, (mpq_ptr)0);
    mpz_clear(f2n);
    return out;
}

/* zeta(-m) = (-1)^m B_{m+1}/(m+1),  for m >= 1 (rational; 0 at even m). */
static Expr* zeta_exact_negint(long m) {
    if (m > ZETA_EXACT_INT_CAP) return NULL;
    if ((m & 1L) == 0) return expr_new_integer(0);     /* even m: B_{odd}=0 */

    mpq_t b, out;
    mpq_inits(b, out, (mpq_ptr)0);
    zbern_get_q(b, (size_t)(m + 1));                   /* B_{m+1} */
    mpq_set_si(out, 1, (unsigned long)(m + 1));         /* 1/(m+1) */
    mpq_mul(out, out, b);
    if (m & 1L) mpq_neg(out, out);                      /* (-1)^m, m odd */
    mpq_canonicalize(out);
    Expr* r = zeta_expr_from_mpq(out);
    mpq_clears(b, out, (mpq_ptr)0);
    return r;
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ */
/* Complex-MPFR toolkit (pairs of mpfr_t; no MPC library available).   */
/* Mirrors the `gcx` helpers in src/gamma.c -- alias-safe, explicit    */
/* working precision `p`.                                              */
/* ------------------------------------------------------------------ */

typedef struct { mpfr_t re, im; } zcx;

static void zcx_init(zcx* z, mpfr_prec_t p) { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void zcx_clear(zcx* z)               { mpfr_clear(z->re);    mpfr_clear(z->im);    }
static void zcx_set(zcx* d, const zcx* s)   { mpfr_set(d->re, s->re, ZRND); mpfr_set(d->im, s->im, ZRND); }

static void zcx_add(zcx* out, const zcx* a, const zcx* b) {
    mpfr_add(out->re, a->re, b->re, ZRND);
    mpfr_add(out->im, a->im, b->im, ZRND);
}

static void zcx_abs(mpfr_t mag, const zcx* z) { mpfr_hypot(mag, z->re, z->im, ZRND); }

static void zcx_mul(zcx* out, const zcx* a, const zcx* b, mpfr_prec_t p) {
    mpfr_t ac, bd, ad, bc;
    mpfr_inits2(p, ac, bd, ad, bc, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, ZRND);
    mpfr_mul(bd, a->im, b->im, ZRND);
    mpfr_mul(ad, a->re, b->im, ZRND);
    mpfr_mul(bc, a->im, b->re, ZRND);
    mpfr_sub(out->re, ac, bd, ZRND);
    mpfr_add(out->im, ad, bc, ZRND);
    mpfr_clears(ac, bd, ad, bc, (mpfr_ptr)0);
}

/* out = a * r, r a real scalar. */
static void zcx_mul_r(zcx* out, const zcx* a, const mpfr_t r) {
    mpfr_mul(out->re, a->re, r, ZRND);
    mpfr_mul(out->im, a->im, r, ZRND);
}

static void zcx_div(zcx* out, const zcx* a, const zcx* b, mpfr_prec_t p) {
    mpfr_t den, t1, t2, nr, ni;
    mpfr_inits2(p, den, t1, t2, nr, ni, (mpfr_ptr)0);
    mpfr_mul(t1, b->re, b->re, ZRND);
    mpfr_mul(t2, b->im, b->im, ZRND);
    mpfr_add(den, t1, t2, ZRND);                 /* |b|^2 */
    mpfr_mul(t1, a->re, b->re, ZRND);
    mpfr_mul(t2, a->im, b->im, ZRND);
    mpfr_add(nr, t1, t2, ZRND);                  /* ac + bd */
    mpfr_mul(t1, a->im, b->re, ZRND);
    mpfr_mul(t2, a->re, b->im, ZRND);
    mpfr_sub(ni, t1, t2, ZRND);                  /* bc - ad */
    mpfr_div(out->re, nr, den, ZRND);
    mpfr_div(out->im, ni, den, ZRND);
    mpfr_clears(den, t1, t2, nr, ni, (mpfr_ptr)0);
}

static void zcx_exp(zcx* out, const zcx* z, mpfr_prec_t p) {
    mpfr_t ea, c, s;
    mpfr_inits2(p, ea, c, s, (mpfr_ptr)0);
    mpfr_exp(ea, z->re, ZRND);
    mpfr_sin_cos(s, c, z->im, ZRND);
    mpfr_mul(out->re, ea, c, ZRND);
    mpfr_mul(out->im, ea, s, ZRND);
    mpfr_clears(ea, c, s, (mpfr_ptr)0);
}

/* out = Log(z), principal branch. */
static void zcx_log(zcx* out, const zcx* z, mpfr_prec_t p) {
    mpfr_t mag, ang;
    mpfr_inits2(p, mag, ang, (mpfr_ptr)0);
    mpfr_hypot(mag, z->re, z->im, ZRND);
    mpfr_atan2(ang, z->im, z->re, ZRND);
    mpfr_log(out->re, mag, ZRND);
    mpfr_set(out->im, ang, ZRND);
    mpfr_clears(mag, ang, (mpfr_ptr)0);
}

/* out = base^expo (principal branch) = exp(expo * Log(base)). */
static void zcx_pow(zcx* out, const zcx* base, const zcx* expo, mpfr_prec_t p) {
    zcx lg, prod;
    zcx_init(&lg, p); zcx_init(&prod, p);
    zcx_log(&lg, base, p);
    zcx_mul(&prod, expo, &lg, p);
    zcx_exp(out, &prod, p);
    zcx_clear(&lg); zcx_clear(&prod);
}

/* ------------------------------------------------------------------ */
/* Euler-Maclaurin Hurwitz zeta kernel                                 */
/* ------------------------------------------------------------------ */

/* Number of head terms N for working precision wp and |s|. The correction
 * series error falls like (2 pi N)^-(2M); taking N ~ digits + |s| keeps M
 * modest and the head sum cheap. */
static unsigned long zeta_em_terms(mpfr_prec_t wp, const zcx* s) {
    double digits = (double)wp * 0.30103 + 5.0;
    double smag = hypot(mpfr_get_d(s->re, ZRND), mpfr_get_d(s->im, ZRND));
    double n = digits + smag + 10.0;
    if (n < 8.0) n = 8.0;
    if (n > 200000.0) n = 200000.0;
    return (unsigned long)n;
}

/* zeta(s, a) for complex s, a (the two-argument Zeta convention) into `out`.
 * Returns false on a non-finite result (e.g. the s = 1 pole).
 *
 * Unlike HurwitzZeta, two-argument Zeta uses the *symmetric* power
 * ((a+k)^2)^(-s/2) (= |a+k|^-s on the real axis), so for Re(a) < 0 the head
 * terms whose base lies in the left half-plane differ from the principal branch
 * (a+k)^-s.  We sum the head with the symmetric power and run the
 * Euler-Maclaurin tail from a+N, where N is forced large enough that Re(a+N) is
 * comfortably positive (there the symmetric and principal powers coincide). */
static bool zeta_hurwitz_cx(zcx* out, const zcx* s, const zcx* a, mpfr_prec_t wp) {
    unsigned long N = zeta_em_terms(wp, s);
    /* Push the tail base a+N safely into the right half-plane for Re(a) < 0. */
    {
        double are = mpfr_get_d(a->re, ZRND);
        if (are < 0) {
            unsigned long need = (unsigned long)(-are) + 8UL;
            if (need > N) N = need;
        }
    }

    zcx S, negs, hnegs, base, base2, term, w, P0, f, winv2, w2, sm1, r, tmp, tmp2;
    zcx_init(&S, wp); zcx_init(&negs, wp); zcx_init(&hnegs, wp);
    zcx_init(&base, wp); zcx_init(&base2, wp); zcx_init(&term, wp);
    zcx_init(&w, wp); zcx_init(&P0, wp); zcx_init(&f, wp); zcx_init(&winv2, wp);
    zcx_init(&w2, wp); zcx_init(&sm1, wp); zcx_init(&r, wp); zcx_init(&tmp, wp);
    zcx_init(&tmp2, wp);

    mpfr_t eps, mag, prevmag, sabs, cj, bq_d;
    mpfr_inits2(wp, eps, mag, prevmag, sabs, cj, bq_d, (mpfr_ptr)0);
    mpfr_set_ui(eps, 1, ZRND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 12 ? wp - 8 : 1), ZRND);

    /* -s and -s/2 */
    mpfr_neg(negs.re, s->re, ZRND);
    mpfr_neg(negs.im, s->im, ZRND);
    mpfr_div_2ui(hnegs.re, negs.re, 1, ZRND);
    mpfr_div_2ui(hnegs.im, negs.im, 1, ZRND);

    /* Head sum S = Sum_{k=0}^{N-1} ((a+k)^2)^(-s/2), skipping any k with
     * a+k == 0 (symmetric power keeps the singular term excluded). */
    mpfr_set_ui(S.re, 0, ZRND); mpfr_set_ui(S.im, 0, ZRND);
    for (unsigned long k = 0; k < N; k++) {
        mpfr_add_ui(base.re, a->re, k, ZRND);
        mpfr_set(base.im, a->im, ZRND);
        zcx_abs(mag, &base);
        if (mpfr_sgn(mag) == 0) continue;          /* k + a == 0 excluded */
        zcx_mul(&base2, &base, &base, wp);         /* (a+k)^2 */
        zcx_pow(&term, &base2, &hnegs, wp);        /* ((a+k)^2)^(-s/2) */
        zcx_add(&S, &S, &term);
    }

    /* w = a + N. */
    mpfr_add_ui(w.re, a->re, N, ZRND);
    mpfr_set(w.im, a->im, ZRND);

    /* P0 = w^-s. */
    zcx_pow(&P0, &w, &negs, wp);

    /* + w^(1-s)/(s-1) = (P0 * w)/(s-1). */
    mpfr_sub_ui(sm1.re, s->re, 1, ZRND);
    mpfr_set(sm1.im, s->im, ZRND);
    zcx_mul(&tmp, &P0, &w, wp);
    zcx_div(&tmp2, &tmp, &sm1, wp);
    zcx_add(&S, &S, &tmp2);

    /* + 1/2 w^-s. */
    mpfr_set_d(cj, 0.5, ZRND);
    zcx_mul_r(&tmp, &P0, cj);
    zcx_add(&S, &S, &tmp);

    /* Correction series.  f_j = w^-(2j-1), starting f_1 = 1/w; advance by
     * winv2 = 1/w^2 each step.  r_j = (s)_{2j-1}, starting r_1 = s. */
    {
        zcx one; zcx_init(&one, wp);
        mpfr_set_ui(one.re, 1, ZRND); mpfr_set_ui(one.im, 0, ZRND);
        zcx_div(&f, &one, &w, wp);                 /* 1/w */
        zcx_mul(&w2, &w, &w, wp);
        zcx_div(&winv2, &one, &w2, wp);            /* 1/w^2 */
        zcx_clear(&one);
    }
    zcx_set(&r, s);                                /* (s)_1 = s */

    mpfr_set_inf(prevmag, 1);
    mpq_t bq;
    mpq_init(bq);
    mpz_t fac;
    mpz_init(fac);

    unsigned long jmax = (unsigned long)wp + 8;
    for (unsigned long j = 1; j <= jmax; j++) {
        /* coeff c_j = B_{2j}/(2j)!. */
        zbern_get_q(bq, (size_t)(2 * j));
        if (mpq_sgn(bq) == 0) {                    /* never (even index), but be safe */
            /* update r, f and continue */
        }
        zeta_factorial(fac, (long)(2 * j));
        mpfr_set_q(cj, bq, ZRND);
        mpfr_set_z(bq_d, fac, ZRND);
        mpfr_div(cj, cj, bq_d, ZRND);              /* B_{2j}/(2j)! */

        /* term = c_j * r * P0 * f. */
        zcx_mul(&tmp, &r, &P0, wp);
        zcx_mul(&tmp2, &tmp, &f, wp);
        zcx_mul_r(&term, &tmp2, cj);
        zcx_abs(mag, &term);

        /* Optimal truncation: stop once the terms start growing again. */
        if (j > 1 && mpfr_cmp(mag, prevmag) > 0) break;
        zcx_add(&S, &S, &term);
        mpfr_set(prevmag, mag, ZRND);

        /* Converged? |term| < eps * |S|. */
        zcx_abs(sabs, &S);
        mpfr_mul(sabs, sabs, eps, ZRND);
        if (mpfr_cmp(mag, sabs) < 0) break;

        /* r *= (s + 2j-1)(s + 2j). */
        mpfr_add_ui(tmp.re, s->re, 2 * j - 1, ZRND); mpfr_set(tmp.im, s->im, ZRND);
        zcx_mul(&tmp2, &r, &tmp, wp);
        mpfr_add_ui(tmp.re, s->re, 2 * j, ZRND);     mpfr_set(tmp.im, s->im, ZRND);
        zcx_mul(&r, &tmp2, &tmp, wp);
        /* f *= 1/w^2. */
        zcx_mul(&tmp, &f, &winv2, wp);
        zcx_set(&f, &tmp);
    }

    mpz_clear(fac);
    mpq_clear(bq);
    zcx_set(out, &S);
    bool ok = mpfr_number_p(out->re) && mpfr_number_p(out->im);

    mpfr_clears(eps, mag, prevmag, sabs, cj, bq_d, (mpfr_ptr)0);
    zcx_clear(&S); zcx_clear(&negs); zcx_clear(&hnegs); zcx_clear(&base);
    zcx_clear(&base2); zcx_clear(&term);
    zcx_clear(&w); zcx_clear(&P0); zcx_clear(&f); zcx_clear(&winv2);
    zcx_clear(&w2); zcx_clear(&sm1); zcx_clear(&r); zcx_clear(&tmp); zcx_clear(&tmp2);
    return ok;
}

/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool zeta_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, ZRND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          ZRND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        ZRND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          ZRND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, ZRND);
        mpfr_div_si(out, out, (long)d, ZRND);
        return true;
    }
    return false;
}

/* Result precision (bits) for Zeta[s, a] under Mathematica contagion: the
 * minimum precision among the inexact leaves of s and a, floored at machine
 * (53). A machine-precision argument forces a machine-precision result even if
 * the other argument is high-precision MPFR. See numeric_min_inexact_bits. */
static mpfr_prec_t zeta_out_prec(const Expr* s, const Expr* a) {
    long ps = numeric_min_inexact_bits(s);
    long pa = numeric_min_inexact_bits(a);
    long m  = (ps && pa) ? (ps < pa ? ps : pa) : (ps ? ps : pa);
    return m < 53 ? 53 : (mpfr_prec_t)m;
}

/* Build the numeric result: Real / Complex[Real,Real] for machine precision,
 * MPFR parts for arbitrary precision.  `real_only` forces a real result. */
static Expr* zeta_result(const mpfr_t re, const mpfr_t im,
                         mpfr_prec_t out_prec, bool real_only) {
    if (out_prec <= 53) {
        Expr* rr = expr_new_real(mpfr_get_d(re, ZRND));
        if (real_only) return rr;
        return make_complex(rr, expr_new_real(mpfr_get_d(im, ZRND)));
    }
    Expr* rr = expr_new_mpfr_bits(out_prec);
    mpfr_set(rr->data.mpfr, re, ZRND);
    if (real_only) return rr;
    Expr* ii = expr_new_mpfr_bits(out_prec);
    mpfr_set(ii->data.mpfr, im, ZRND);
    return make_complex(rr, ii);
}

/* Decompose a numeric value into borrowed (re, im). A plain real yields
 * (value, shared zero); a Complex[..] yields its two parts. */
static void zeta_decompose(Expr* v, Expr* zero, Expr** re, Expr** im) {
    Expr *r, *i;
    if (is_complex(v, &r, &i)) { *re = r; *im = i; }
    else { *re = v; *im = zero; }
}

/* Numeric Zeta[s, a] (a defaults to 1 for the Riemann case). At least one of
 * s, a must carry an inexact part. Returns Real/Complex/MPFR, ComplexInfinity
 * at the pole, or NULL if the parts are not numeric. */
static Expr* zeta_numeric(Expr* s, Expr* a) {
    Expr* zero = expr_new_integer(0);
    Expr *sre, *sim, *are, *aim;
    zeta_decompose(s, zero, &sre, &sim);
    zeta_decompose(a, zero, &are, &aim);

    mpfr_prec_t out_prec = zeta_out_prec(s, a);
    mpfr_prec_t wp = out_prec + 64;

    /* Guard bits for Re(s) < 0 (mirrors hz_numeric in hurwitzzeta.c): the
     * head terms grow like (N+|a|)^(-Re s) while the value is far smaller, so
     * the Euler-Maclaurin sum cancels ~(-Re s) log2(N+|a|) bits.  Without them
     * Zeta[-38.5, 2.3] came out as -1.45e22 (true value -1.72e14) at every
     * precision.  `nest` bounds the N zeta_em_terms picks once the guard is
     * added (it also covers the Re(a) < 0 extension of N). */
    {
        double sre_d = 0.0, are_d = 0.0, aim_d = 0.0;
        mpfr_t t;
        mpfr_init2(t, 64);
        if (zeta_set_mpfr(t, sre)) sre_d = mpfr_get_d(t, ZRND);
        if (zeta_set_mpfr(t, are)) are_d = mpfr_get_d(t, ZRND);
        if (zeta_set_mpfr(t, aim)) aim_d = mpfr_get_d(t, ZRND);
        mpfr_clear(t);
        if (sre_d < 0.0 && isfinite(sre_d)) {
            double digits = (double)wp * 0.30103 + 5.0;
            double nest = digits + 5.0 * fabs(sre_d) + 74.0;
            double guard = -sre_d * log2(nest + hypot(are_d, aim_d) + 9.0) + 16.0;
            if (guard < 1e7) wp += (mpfr_prec_t)guard;
        }
    }

    bool real_only = (sim == zero) && (aim == zero);

    zcx sc, ac, res;
    zcx_init(&sc, wp); zcx_init(&ac, wp); zcx_init(&res, wp);
    bool ok = zeta_set_mpfr(sc.re, sre) && zeta_set_mpfr(sc.im, sim) &&
              zeta_set_mpfr(ac.re, are) && zeta_set_mpfr(ac.im, aim);
    Expr* out = NULL;
    if (ok) {
        if (zeta_hurwitz_cx(&res, &sc, &ac, wp)) {
            out = zeta_result(res.re, res.im, out_prec, real_only);
        } else {
            out = expr_new_symbol(SYM_ComplexInfinity);
        }
    }
    zcx_clear(&sc); zcx_clear(&ac); zcx_clear(&res);
    expr_free(zero);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Zeta[s]  (Riemann zeta)                                             */
/* ------------------------------------------------------------------ */

static Expr* zeta_one_arg(Expr* s) {
    /* 1. Exact integer arguments -> closed forms. */
    long n;
    if (zeta_exact_int(s, &n)) {
        if (n == 1) return expr_new_symbol(SYM_ComplexInfinity);
        if (n == 0) return make_rational(-1, 2);
        if (n > 0) {
            if ((n & 1L) == 0) {                    /* even positive: rational*Pi^n */
                Expr* ev = zeta_exact_even(n);
                if (ev) return ev;
            }
            /* odd positive (>=3): no closed form -> symbolic (unless numeric) */
        } else {                                    /* n < 0 */
            Expr* nv = zeta_exact_negint(-n);
            if (nv) return nv;
        }
    }

    /* 2. Symbolic limiting value. */
    if (zeta_is_symbol(s, "Infinity")) return expr_new_integer(1);

#ifdef USE_MPFR
    /* 2.5. Arbitrary-precision real s, and complex s, via arb (acb_dirichlet_zeta):
     * far faster than mpfr_zeta and the hand-rolled complex Euler-Maclaurin kernel
     * below (e.g. ~200x at 50 digits). Machine-double real s is left to the direct
     * mpfr_zeta -> double path (one library call, no arb ball overhead). A decline
     * (NULL at the s=1 pole) falls through to the paths below, which yield
     * ComplexInfinity. */
    {
        Expr *cre, *cim;
        /* Only ARBITRARY-precision (MPFR) real or complex s routes to arb; a
         * machine Real (or machine complex with Real parts) keeps the existing
         * mpfr_zeta -> double / zeta_numeric -> double paths so machine-in stays
         * machine-out (the arb bridge returns EXPR_MPFR). */
        bool cplx_mpfr = is_complex(s, &cre, &cim) &&
                         ((cre && cre->type == EXPR_MPFR) ||
                          (cim && cim->type == EXPR_MPFR));
        if (s->type == EXPR_MPFR || cplx_mpfr) {
            Expr* v = flint_num_zeta(s);
            if (v) return v;
        }
    }

    /* 3/4. Numeric real Riemann zeta via MPFR (machine or arbitrary). */
    if (s->type == EXPR_REAL || s->type == EXPR_MPFR) {
        mpfr_prec_t prec = (s->type == EXPR_MPFR) ? mpfr_get_prec(s->data.mpfr) : 53;
        mpfr_t sv, rv;
        mpfr_init2(sv, prec + 16); mpfr_init2(rv, prec + 16);
        zeta_set_mpfr(sv, s);
        mpfr_zeta(rv, sv, ZRND);
        Expr* out;
        if (!mpfr_number_p(rv)) {
            out = expr_new_symbol(SYM_ComplexInfinity);   /* pole at s = 1 */
        } else if (s->type == EXPR_MPFR) {
            out = expr_new_mpfr_bits(prec);
            mpfr_set(out->data.mpfr, rv, ZRND);
        } else {
            out = expr_new_real(mpfr_get_d(rv, ZRND));
        }
        mpfr_clear(sv); mpfr_clear(rv);
        return out;
    }

    /* 5. Complex argument (Complex[..] with an inexact part) -> EM kernel
     *    with a = 1. */
    {
        Expr *re, *im;
        if (is_complex(s, &re, &im) &&
            (zeta_is_inexact(re) || zeta_is_inexact(im))) {
            Expr* one = expr_new_integer(1);
            Expr* out = zeta_numeric(s, one);
            expr_free(one);
            if (out) return out;
        }
    }
#endif

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Zeta[s, a]  (Hurwitz zeta)                                          */
/* ------------------------------------------------------------------ */

/* Zeta[s, 1/2] = (2^s - 1) Zeta[s].  Valid because a = 1/2 > 0, so the symmetric
 * power ((k+1/2)^2)^(-s/2) coincides with the principal (k+1/2)^-s. */
static Expr* zeta_half(Expr* s) {
    Expr* twos = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_new_integer(2), expr_copy(s) }, 2);
    Expr* coeff = expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ twos, expr_new_integer(-1) }, 2);
    Expr* zs = expr_new_function(expr_new_symbol(SYM_Zeta),
                     (Expr*[]){ expr_copy(s) }, 1);
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ coeff, zs }, 2);
    return eval_and_free(prod);
}

/* out = v, for a non-negative int64 v (unsigned long may be 32 bits). */
static void zeta_mpz_set_u64(mpz_t out, uint64_t v) {
    mpz_import(out, 1, 1, sizeof v, 0, 0, &v);
}

/* acc += sign (p/q)^-s for p, q > 0 and an exact integer s: the head terms of
 * a shifted Hurwitz sum, kept as one exact rational so that Zeta[2, -200] is a
 * GMP sum and a two-term Plus rather than 200 evaluated powers. */
static void zeta_head_add(mpq_t acc, int64_t p, int64_t q, long s, int sign) {
    unsigned long e = (unsigned long)(s >= 0 ? s : -s);
    mpz_t num, den;
    mpz_inits(num, den, (mpz_ptr)0);
    zeta_mpz_set_u64(s >= 0 ? den : num, (uint64_t)p);
    zeta_mpz_set_u64(s >= 0 ? num : den, (uint64_t)q);
    mpz_pow_ui(num, num, e);
    mpz_pow_ui(den, den, e);
    mpq_t t;
    mpq_init(t);
    mpq_set_num(t, num);
    mpq_set_den(t, den);
    mpq_canonicalize(t);
    if (sign < 0) mpq_sub(acc, acc, t); else mpq_add(acc, acc, t);
    mpq_clear(t);
    mpz_clears(num, den, (mpz_ptr)0);
}

/* rest + h, evaluated; `rest` is consumed. */
static Expr* zeta_plus_rational(Expr* rest, const mpq_t h) {
    mpz_t num, den;
    mpz_inits(num, den, (mpz_ptr)0);
    mpq_get_num(num, h);
    mpq_get_den(den, h);
    Expr* r = make_rational_mpz(num, den);
    mpz_clears(num, den, (mpz_ptr)0);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                                           (Expr*[]){ r, rest }, 2));
}

/* True when Zeta[s] has a closed form: s an exact integer other than an odd
 * s >= 3. Integer-a shifts past the trivial a in {-1, 0, 1, 2} (whose extra
 * terms are all 1^-s) are taken only then, so Zeta[3, 4] stays folded. */
static bool zeta_closed_int_s(Expr* s) {
    long n;
    return zeta_exact_int(s, &n) && !(n >= 3 && (n & 1L));
}

/* Exact Hurwitz at a positive integer a = m:
 *   zeta(s, m) = Zeta[s] - Sum_{k=1}^{m-1} k^-s.
 * For any s when m <= 2 (Zeta[x, 2] -> -1 + Zeta[x]); past that only when
 * Zeta[s] has a closed form (Zeta[4,5] -> Pi^4/90 - 22369/20736, while
 * Zeta[3, 4] and Zeta[x, 3] stay). Returns NULL past the cap. */
static Expr* zeta_hurwitz_int_a(Expr* s, long m) {
    if (m < 1 || m > ZETA_HURWITZ_A_CAP) return NULL;
    if (m >= 3 && !zeta_closed_int_s(s)) return NULL;

    Expr* zs = expr_new_function(expr_new_symbol(SYM_Zeta), (Expr*[]){ expr_copy(s) }, 1);
    if (m == 1) return eval_and_free(zs);          /* zeta(s, 1) = zeta(s) */
    long si;
    if (zeta_exact_int(s, &si)) {                  /* exact rational head */
        mpq_t h;
        mpq_init(h);
        for (long k = 1; k < m; k++) zeta_head_add(h, k, 1, si, -1);
        Expr* out = zeta_plus_rational(eval_and_free(zs), h);
        mpq_clear(h);
        return out;
    }

    /* Build Zeta[s] - Sum_{k=1}^{m-1} k^-s. */
    size_t terms = (size_t)m;                       /* zs + (m-1) subtracted terms */
    Expr** parts = (Expr**)malloc(terms * sizeof(Expr*));
    if (!parts) { expr_free(zs); return NULL; }
    parts[0] = zs;
    for (long k = 1; k < m; k++) {
        /* -(k^-s) = Times[-1, Power[k, -s]]. */
        Expr* negs = expr_new_function(expr_new_symbol(SYM_Times),
                         (Expr*[]){ expr_new_integer(-1), expr_copy(s) }, 2);
        Expr* kpow = expr_new_function(expr_new_symbol(SYM_Power),
                         (Expr*[]){ expr_new_integer(k), negs }, 2);
        parts[(size_t)k] = expr_new_function(expr_new_symbol(SYM_Times),
                         (Expr*[]){ expr_new_integer(-1), kpow }, 2);
    }
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), parts, terms);
    free(parts);
    return eval_and_free(sum);
}

/* Plus[parts[0..n-1]] evaluated; `parts` and its entries are consumed. */
static Expr* zeta_plus_free(Expr** parts, size_t n) {
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), parts, n);
    free(parts);
    return eval_and_free(sum);
}

/* b^-s for an exact rational b > 0, unevaluated (Power[b, Times[-1, s]]). */
static Expr* zeta_pow_negs(Expr* b, Expr* s) {
    Expr* negs = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ expr_new_integer(-1), expr_copy(s) }, 2);
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ b, negs }, 2);
}

/* Zeta[s, r] for the shift target r, evaluated; NULL (freed) when it stays
 * unevaluated, which is when a shifted form would only restate the input. */
static Expr* zeta_shift_target(Expr* s, Expr* r) {
    Expr* v = eval_and_free(expr_new_function(expr_new_symbol(SYM_Zeta),
                  (Expr*[]){ expr_copy(s), r }, 2));
    if (v->type == EXPR_FUNCTION && v->data.function.head->type == EXPR_SYMBOL &&
        v->data.function.head->data.symbol.name == SYM_Zeta) {
        expr_free(v);
        return NULL;
    }
    return v;
}

/* Exact Hurwitz at a rational a > 1, a not an integer, for an exact integer s:
 * shift down to r = a - m in (0, 1),
 *   zeta(s, a) = zeta(s, r) - Sum_{k<m} (r + k)^-s,
 * when Zeta[s, r] evaluates (r = 1/2). NULL otherwise. */
static Expr* zeta_hurwitz_shift_down(Expr* s, long si, int64_t n, int64_t d) {
    int64_t m = n / d;
    if (m > ZETA_HURWITZ_A_CAP) return NULL;
    int64_t rn = n - m * d;
    Expr* shifted = zeta_shift_target(s, make_rational(rn, d));
    if (!shifted) return NULL;
    mpq_t h;
    mpq_init(h);
    for (int64_t k = 0; k < m; k++) zeta_head_add(h, rn + k * d, d, si, -1);
    Expr* out = zeta_plus_rational(shifted, h);
    mpq_clear(h);
    return out;
}

/* Exact Hurwitz at a <= 0 (integer or rational), symmetric convention:
 *   zeta(s, a) = Sum_{0<=k<m, k+a<0} |k+a|^-s  [+ 1 when s = 0, a integer]
 *                + zeta(s, a + m),
 * with m the number of k for which k + a <= 0. For integer a = -m the shifted
 * term is Zeta[s] (the dropped k = m term is 0^-s); otherwise a + m = r lies in
 * (0, 1). Expanded where the reference expands: a = 0 or -1 for any s; a <= -2
 * when Zeta[s] has a closed form; non-integer a for an exact integer s whose shifted
 * Zeta[s, r] evaluates. NULL otherwise, or past ZETA_HURWITZ_A_CAP terms.
 * A rational a > 1 is shifted down the same way (zeta_hurwitz_shift_down). */
static Expr* zeta_hurwitz_shift(Expr* s, Expr* a) {
    int64_t n, d;
    if (!is_rational(a, &n, &d)) return NULL;
    long sint;
    bool s_int = zeta_exact_int(s, &sint);
    if (n > 0) {
        if (d > 1 && n > d && s_int) return zeta_hurwitz_shift_down(s, sint, n, d);
        return NULL;
    }

    if (d == 1) {                                   /* a = -m */
        int64_t m = -n;
        if (m >= 2 && !zeta_closed_int_s(s)) return NULL;
        if (m > ZETA_HURWITZ_A_CAP) return NULL;
        bool s_zero = s_int && sint == 0;
        if (s_int) {                                /* exact rational head */
            mpq_t h;
            mpq_init(h);
            for (int64_t j = 1; j <= m; j++) zeta_head_add(h, j, 1, sint, 1);
            if (s_zero) mpq_set_ui(h, (unsigned long)m + 1, 1);   /* m ones, and 0^0 = 1 */
            Expr* out = zeta_plus_rational(eval_and_free(expr_new_function(
                            expr_new_symbol(SYM_Zeta), (Expr*[]){ expr_copy(s) }, 1)), h);
            mpq_clear(h);
            return out;
        }
        size_t terms = (size_t)m + 1 + (s_zero ? 1 : 0);
        Expr** parts = (Expr**)malloc(terms * sizeof(Expr*));
        if (!parts) return NULL;
        size_t t = 0;
        parts[t++] = expr_new_function(expr_new_symbol(SYM_Zeta),
                                       (Expr*[]){ expr_copy(s) }, 1);
        for (int64_t j = 1; j <= m; j++)
            parts[t++] = zeta_pow_negs(expr_new_integer(j), s);
        if (s_zero) parts[t++] = expr_new_integer(1);   /* 0^0 */
        return zeta_plus_free(parts, t);
    }

    /* a = n/d, d > 1: m = ceil(-a) terms with k + a < 0, then r = a + m. */
    if (!s_int) return NULL;
    int64_t m = (-n + d - 1) / d;
    if (m > ZETA_HURWITZ_A_CAP) return NULL;
    Expr* shifted = zeta_shift_target(s, make_rational(n + m * d, d));
    if (!shifted) return NULL;
    mpq_t h;
    mpq_init(h);
    for (int64_t k = 0; k < m; k++)                 /* |k + a| = (-n - k d)/d */
        zeta_head_add(h, -n - k * d, d, sint, 1);
    Expr* out = zeta_plus_rational(shifted, h);
    mpq_clear(h);
    return out;
}

static Expr* zeta_two_arg(Expr* s, Expr* a) {
    /* 1. zeta(s, 1) = zeta(s). */
    if (a->type == EXPR_INTEGER && a->data.integer == 1)
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Zeta),
                             (Expr*[]){ expr_copy(s) }, 1));

    /* 2. a == 1/2: Zeta[s, 1/2] = (2^s - 1) Zeta[s] (exact s only). */
    if (!zeta_is_inexact(s) && zeta_is_half(a)) return zeta_half(s);

    /* 3. Exact Hurwitz at a positive integer a (only when s is not inexact;
     *    inexact s falls through to the numeric kernel). */
    if (!zeta_is_inexact(s)) {
        long m;
        if (zeta_exact_int(a, &m) && m >= 1) {
            Expr* closed = zeta_hurwitz_int_a(s, m);
            if (closed) return closed;
        }
        /* 3a. a <= 0, or a rational a > 1, shifted into (0, 1]. */
        Expr* np = zeta_hurwitz_shift(s, a);
        if (np) return np;
        /* 3b. zeta(0, a) = 1/2 - a for a rational a > 0. */
        int64_t n, d;
        if (s->type == EXPR_INTEGER && s->data.integer == 0 &&
            is_rational(a, &n, &d) && n > 0)
            return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                       (Expr*[]){ make_rational(1, 2),
                                  expr_new_function(expr_new_symbol(SYM_Times),
                                      (Expr*[]){ expr_new_integer(-1), expr_copy(a) }, 2) }, 2));
    }

#ifdef USE_MPFR
    /* 4. Numeric (real or complex) Hurwitz: at least one inexact operand. */
    {
        Expr* zero = expr_new_integer(0);
        Expr *sre, *sim, *are, *aim;
        zeta_decompose(s, zero, &sre, &sim);
        zeta_decompose(a, zero, &are, &aim);
        bool any_inexact = zeta_is_inexact(sre) || zeta_is_inexact(sim) ||
                           zeta_is_inexact(are) || zeta_is_inexact(aim);
        expr_free(zero);
        if (any_inexact) {
            Expr* out = zeta_numeric(s, a);
            if (out) return out;
        }
    }
#endif

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

static Expr* zeta_emit_argt(size_t argc) {
    mth_message("Zeta", "argt", "Zeta called with %zu argument%s; 1 or 2 arguments are expected.", argc, argc == 1 ? "" : "s");
    return NULL;
}

/* True if `e` is a two-argument Rule or RuleDelayed (a trailing option). */
static bool zeta_is_rule(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
           e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rule ||
            e->data.function.head->data.symbol.name == SYM_RuleDelayed);
}

/* The real part of a numeric `e` as a double, and whether `e` is exact. False
 * when `e` is not a number (or a complex number with a non-numeric part). */
static bool zeta_re_part(const Expr* e, double* re, bool* exact) {
    Expr *r, *i;
    const Expr* x = is_complex((Expr*)e, &r, &i) ? r : e;
    *exact = !zeta_is_inexact(x) && !(x != e && zeta_is_inexact(i));
#ifdef USE_MPFR
    mpfr_t t;
    mpfr_init2(t, 64);
    bool ok = zeta_set_mpfr(t, x);
    if (ok) *re = mpfr_get_d(t, ZRND);
    mpfr_clear(t);
    return ok;
#else
    if (x->type == EXPR_REAL) { *re = x->data.real; return true; }
    int64_t n, d;
    if (is_rational(x, &n, &d)) { *re = (double)n / (double)d; return true; }
    return false;
#endif
}

/* Zeta[s, a, IncludeSingularTerm -> True] at an integer-valued a <= 0: the
 * k = -a term is 0^-s. Sets *out to ComplexInfinity when Re s > 0 and returns
 * true; returns true with *out NULL for an exact non-integer s with Re s < 0
 * (stays unevaluated) or a non-numeric s; false when the option changes
 * nothing (non-integer a, s = 0, a negative integer s, or inexact Re s <= 0). */
static bool zeta_singular_term(Expr* s, Expr* a, Expr** out) {
    *out = NULL;
    double av, sv;
    bool aex, sex;
    Expr *ar, *ai;
    if (is_complex(a, &ar, &ai) || !zeta_re_part(a, &av, &aex)) return false;
    if (!(av <= 0.0) || av != floor(av)) return false;
    if (!zeta_re_part(s, &sv, &sex)) return true;          /* symbolic s */
    if (sv > 0.0) { *out = expr_new_symbol(SYM_ComplexInfinity); return true; }
    long si;
    if (sex && !zeta_exact_int(s, &si)) return true;       /* exact, non-integer */
    return false;
}

Expr* builtin_zeta(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    /* Trailing options: validate, apply IncludeSingularTerm, then evaluate the
     * option-free call (so every numeric path below sees plain arguments). */
    if (argc > 2 || (argc == 2 && zeta_is_rule(args[1]))) {
        const char* names[1] = { SYM_IncludeSingularTerm };
        bool ist = false;
        size_t npos;
        if (!options_parse_bools(res, "Zeta", 2, names, &ist, 1, &npos)) return NULL;
        if (npos == 0) return zeta_emit_argt(argc);
        if (ist && npos == 2) {
            /* A buffer argument: thread element by element, so each integer
             * a <= 0 meets the option (the kernel knows nothing of it). */
            if (is_ndarray(args[0]) || is_ndarray(args[1]))
                return ndarray_delist_and_reeval(res);
            Expr* sing;
            if (zeta_singular_term(args[0], args[1], &sing)) return sing;
        }
        Expr* plain = expr_new_function(expr_new_symbol(SYM_Zeta),
                          (Expr*[]){ expr_copy(args[0]), npos == 2 ? expr_copy(args[1]) : NULL },
                          npos);
        Expr* v = eval_and_free(plain);
        if (v->type == EXPR_FUNCTION && v->data.function.head->type == EXPR_SYMBOL &&
            v->data.function.head->data.symbol.name == SYM_Zeta) {
            expr_free(v);                                  /* stays unevaluated */
            return NULL;
        }
        return v;
    }

    if (argc == 1 && is_interval(args[0])) {
        Expr* r = interval_apply_function("Zeta", args[0]); if (r) return r;
    }
    /* Machine-precision call: the raw-double kernels shared with the NDArray /
     * packed / Compile[] surfaces (sf_machine_zeta.c), so all surfaces agree
     * bit-for-bit. For a real a <= 0 the kernel shifts a into (0, 1] in the
     * symmetric convention. A decline (pole, domain, overflow) falls through
     * to the arbitrary-precision paths unchanged. */
    {
        SFMachineArgs ma;
        if ((argc == 1 || argc == 2) && sf_machine_scalar_args(res, &ma)) {
            double vr, vi = 0.0;
            bool ok;
            if (argc == 1)
                ok = ma.any_complex ? sf_zeta_c(ma.re[0], ma.im[0], &vr, &vi)
                                    : sf_zeta_r(ma.re[0], &vr);
            else
                ok = ma.any_complex
                         ? sf_zeta2_c(ma.re[0], ma.im[0], ma.re[1], ma.im[1], &vr, &vi)
                         : sf_zeta2_r(ma.re[0], ma.re[1], &vr);
            if (ok) return ma.any_complex ? sf_machine_complex_result(vr, vi)
                                          : expr_new_real(vr);
        }
    }
    if (argc == 1) return zeta_one_arg(args[0]);
    if (argc == 2) return zeta_two_arg(args[0], args[1]);
    return zeta_emit_argt(argc);                           /* argc == 0 */
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

/* Derivative[n][Zeta][s] at an inexact real or complex s: zeta^(n)(s) from
 * arb's Taylor jet (flint_num_zeta_deriv).  D[Zeta[s], s] stays the inert
 * Derivative[1][Zeta][s], as in the reference; this gives it a value, so
 * -Zeta'[s]/Zeta[s] (Sum of MangoldtLambda[k]/k^s), Exp[-Zeta'[2]] and a
 * Newton step on Zeta are numbers.  A machine argument answers a machine
 * number (the jet is computed at 53 bits plus guard bits).  Registered with
 * the derivative registry (deriv.h): a system rule, not a SubValue. */
static Expr* zeta_derivative_value(const Expr* call) {
    const Expr* d = call->data.function.head->data.function.head;   /* Derivative[n] */
    if (d->data.function.arg_count != 1 || d->data.function.args[0]->type != EXPR_INTEGER
        || call->data.function.arg_count != 1)
        return NULL;
    int64_t n = d->data.function.args[0]->data.integer;
    Expr* s = call->data.function.args[0];
    Expr *cre, *cim;
    bool num = s->type == EXPR_REAL
#ifdef USE_MPFR
               || s->type == EXPR_MPFR
#endif
               || (is_complex(s, &cre, &cim) && arg_is_inexact(s));
    if (!num || n < 1 || n > 64) return NULL;
    Expr* v = flint_num_zeta_deriv((long)n, s);
    if (!v) return NULL;
    if (numeric_min_inexact_bits(s) <= 53) {              /* machine in, machine out */
        Expr* m = numericalize(v, numeric_machine_spec());
        expr_free(v);
        v = m;
    }
    return v;
}

void zeta_init(void) {
    deriv_register_numeric(SYM_Zeta, zeta_derivative_value, 1);
    symtab_add_builtin("Zeta", builtin_zeta);
    symtab_get_def("Zeta")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_options("Zeta", expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_new_symbol(SYM_IncludeSingularTerm),
                       expr_new_symbol(SYM_False) }, 2) }, 1));
    /* Docstring lives in info.c (info_init). */
}

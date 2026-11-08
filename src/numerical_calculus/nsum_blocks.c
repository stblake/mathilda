/*
 * nsum_blocks.c — summation of an integer-indexed series from its terms alone
 *
 * See nsum_blocks.h for the model.  The pieces, bottom up:
 *
 *   nsb_term       a term of a level: the summand itself (level 0) or one
 *                  block sum of the level below (level >= 1).
 *   nsb_interp     the smooth continuation of a level's terms inside one block,
 *                  by a local Lagrange interpolant whose width is chosen per
 *                  point from the (x-scale) smoothness bound.
 *   nsb_block      Sum_{n=A}^{B-1} of a level whose terms are smooth on [A, B):
 *                  short ranges term by term, long ones by Euler–Maclaurin with
 *                  the integral in log n (Gauss–Legendre, doubling until two
 *                  rules agree) and the endpoint derivatives from finite
 *                  differences (Stirling-number expansion of D = log(1 + Delta)).
 *   nsb_classify   the head of a level: SMOOTH, BLOCKED (base b, offset c,
 *                  verified on every breakpoint the head can see), or ROUGH.
 *   nsb_sum_level  one level: Levin for SMOOTH, the head + the next level for
 *                  BLOCKED, with every unsummed remainder bounded honestly.
 */
#include <stdio.h>          /* before mpfr.h: mpfr_fprintf */
#include "nsum_blocks.h"

#ifdef USE_MPFR

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "mtime.h"
#include "seqaccel.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NSB_MAX_DEPTH   4      /* block levels below the summand                 */
#define NSB_QMAX        64     /* widest interpolation stencil                   */
#define NSB_QE_MAX      64     /* widest endpoint finite-difference stencil      */
#define NSB_JMAX        28     /* Euler–Maclaurin correction terms               */
#define NSB_H0          1100   /* level-0 head: three decades of breakpoints     */
#define NSB_H0_FIRST    128    /* ... looked at first                             */
#define NSB_CHEAP_TERM_S 2e-5  /* a term this cheap (s) widens a smooth head      */
#define NSB_DIRECT_MAX  160    /* a range this short is summed term by term      */
#define NSB_GL_RULES    5      /* Gauss–Legendre rules 4, 8, 16, 32, 64          */
#define NSB_FLAG_BIG    0.05   /* log-jump that marks a breakpoint               */
#define NSB_FLAG_SMALL  0.02   /* ... while each side stays on its smooth law    */

typedef struct NsbLevel {
    int                     depth;
    const struct NsbLevel*  parent;   /* level whose blocks this level's terms sum */
    mpz_t                   b;        /* the parent's block base ...               */
    long                    c;        /* ... and offset: term k = [b^(k-1)+c, b^k+c) */
    double                  lmax;     /* log2 of the largest index the budget reaches */
    long                    tq;       /* bits this level's terms must carry          */
} NsbLevel;

typedef struct {
    NsbTermFn    f;
    void*        ud;
    long         tb;              /* target bits                                   */
    long         tq;              /* bits each term must carry: tb plus the
                                   * amplification of the Levin transform over
                                   * the terms (~2^L for L terms of a
                                   * logarithmically convergent series)           */
    mpfr_prec_t  p;               /* working bits                                  */
    double       t_end;           /* budget deadline (monotone seconds)            */
    size_t       max_index_bits;  /* largest index the summand is asked for        */
    bool         budget_hit;
    bool         bad;             /* a non-real / non-numeric term: decline        */
    int          max_level;       /* deepest level reached                         */
    long         gq;              /* accuracy goal (bits) of the term in progress  */
    long         nevals;          /* summand evaluations (diagnostics)             */
    /* Gauss–Legendre nodes/weights on [-1, 1], built on first use. */
    int          gl_n[NSB_GL_RULES];
    mpfr_t*      gl_x[NSB_GL_RULES];
    mpfr_t*      gl_w[NSB_GL_RULES];
    /* B_{2j}/(2j)!, j = 1..NSB_JMAX. */
    mpfr_t       bern[NSB_JMAX + 1];
    /* Signed Stirling numbers of the first kind s(m, r), m, r < NSB_QE_MAX. */
    mpz_t*       stir;
} NsbCtx;

/* ------------------------------------------------------------------ *
 *  Small helpers                                                      *
 * ------------------------------------------------------------------ */

/* log2 |v| as a double (-INFINITY for zero); exact in range for any exponent. */
static double nsb_log2(const mpfr_t v) {
    if (mpfr_zero_p(v)) return -INFINITY;
    long e;
    double d = mpfr_get_d_2exp(&e, v, MPFR_RNDN);
    return (double)e + log2(fabs(d));
}

static double nsb_log2_z(const mpz_t z) {
    if (mpz_sgn(z) == 0) return -INFINITY;
    long e;
    double d = mpz_get_d_2exp(&e, z);
    return (double)e + log2(fabs(d));
}

/* |v| as a double, saturating instead of underflowing to a misleading 0 only
 * when v is genuinely tiny (callers use it for error bookkeeping). */
static double nsb_absd(const mpfr_t v) {
    return fabs(mpfr_get_d(v, MPFR_RNDU));
}

/* MATHILDA_NSB_DEBUG=1 traces each level's shape and estimate on stderr. */
static int nsb_debug(void) {
    static int on = -1;
    if (on < 0) { const char* e = getenv("MATHILDA_NSB_DEBUG"); on = (e && *e && *e != '0') ? atoi(e) : 0; }
    return on;
}

static bool nsb_over_budget(NsbCtx* x) {
    if (x->budget_hit) return true;
    if (mth_mono_seconds() > x->t_end) { x->budget_hit = true; return true; }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Tables: Gauss–Legendre, Bernoulli, Stirling                        *
 * ------------------------------------------------------------------ */

/* Nodes and weights of the n-point Gauss–Legendre rule by Newton's method on
 * P_n, from the classical cosine guess, at the working precision. */
static bool nsb_gl_build(NsbCtx* x, int r) {
    if (x->gl_x[r]) return true;
    int n = x->gl_n[r];
    mpfr_t* X = malloc(sizeof(mpfr_t) * (size_t)n);
    mpfr_t* W = malloc(sizeof(mpfr_t) * (size_t)n);
    if (!X || !W) { free(X); free(W); return false; }
    mpfr_prec_t p = x->p + 16;
    mpfr_t z, p0, p1, p2, dp, t, u;
    mpfr_inits2(p, z, p0, p1, p2, dp, t, u, (mpfr_ptr)0);
    for (int i = 0; i < n; i++) {
        mpfr_init2(X[i], x->p); mpfr_init2(W[i], x->p);
        mpfr_set_d(z, cos(M_PI * (i + 0.75) / (n + 0.5)), MPFR_RNDN);
        for (int it = 0; it < 64; it++) {
            /* P_n(z) by the three-term recurrence; P_n' from P_n, P_{n-1}. */
            mpfr_set_ui(p0, 1, MPFR_RNDN);
            mpfr_set(p1, z, MPFR_RNDN);
            for (int k = 2; k <= n; k++) {
                /* p2 = ((2k-1) z p1 - (k-1) p0) / k */
                mpfr_mul(t, z, p1, MPFR_RNDN);
                mpfr_mul_ui(t, t, (unsigned long)(2 * k - 1), MPFR_RNDN);
                mpfr_mul_ui(u, p0, (unsigned long)(k - 1), MPFR_RNDN);
                mpfr_sub(t, t, u, MPFR_RNDN);
                mpfr_div_ui(p2, t, (unsigned long)k, MPFR_RNDN);
                mpfr_set(p0, p1, MPFR_RNDN);
                mpfr_set(p1, p2, MPFR_RNDN);
            }
            /* dp = n (z p1 - p0) / (z^2 - 1) */
            mpfr_mul(t, z, p1, MPFR_RNDN);
            mpfr_sub(t, t, p0, MPFR_RNDN);
            mpfr_mul_ui(t, t, (unsigned long)n, MPFR_RNDN);
            mpfr_sqr(u, z, MPFR_RNDN);
            mpfr_sub_ui(u, u, 1, MPFR_RNDN);
            mpfr_div(dp, t, u, MPFR_RNDN);
            mpfr_div(t, p1, dp, MPFR_RNDN);
            mpfr_sub(z, z, t, MPFR_RNDN);
            if (mpfr_zero_p(t) || mpfr_get_exp(t) < -(mpfr_exp_t)(p + 4)) break;
        }
        mpfr_set(X[i], z, MPFR_RNDN);
        /* w = 2 / ((1 - z^2) dp^2) */
        mpfr_sqr(u, z, MPFR_RNDN);
        mpfr_ui_sub(u, 1, u, MPFR_RNDN);
        mpfr_sqr(t, dp, MPFR_RNDN);
        mpfr_mul(u, u, t, MPFR_RNDN);
        mpfr_ui_div(W[i], 2, u, MPFR_RNDN);
    }
    mpfr_clears(z, p0, p1, p2, dp, t, u, (mpfr_ptr)0);
    x->gl_x[r] = X; x->gl_w[r] = W;
    return true;
}

/* B_{2j}/(2j)! = (-1)^(j+1) 2 zeta(2j) / (2 pi)^(2j). */
static void nsb_bern_build(NsbCtx* x) {
    mpfr_t tp, z;
    mpfr_inits2(x->p + 16, tp, z, (mpfr_ptr)0);
    mpfr_const_pi(tp, MPFR_RNDN);
    mpfr_mul_ui(tp, tp, 2, MPFR_RNDN);
    for (int j = 1; j <= NSB_JMAX; j++) {
        mpfr_init2(x->bern[j], x->p);
        mpfr_zeta_ui(z, (unsigned long)(2 * j), MPFR_RNDN);
        mpfr_mul_ui(z, z, 2, MPFR_RNDN);
        mpfr_t pw; mpfr_init2(pw, x->p + 16);
        mpfr_pow_ui(pw, tp, (unsigned long)(2 * j), MPFR_RNDN);
        mpfr_div(x->bern[j], z, pw, MPFR_RNDN);
        mpfr_clear(pw);
        if (!(j & 1)) mpfr_neg(x->bern[j], x->bern[j], MPFR_RNDN);
    }
    mpfr_init2(x->bern[0], 2);
    mpfr_clears(tp, z, (mpfr_ptr)0);
}

/* s(m, r) by s(m+1, r) = s(m, r-1) - m s(m, r). */
static bool nsb_stir_build(NsbCtx* x) {
    size_t N = NSB_QE_MAX;
    x->stir = malloc(sizeof(mpz_t) * N * N);
    if (!x->stir) return false;
    for (size_t i = 0; i < N * N; i++) mpz_init(x->stir[i]);
    mpz_set_ui(x->stir[0], 1);
    for (size_t m = 0; m + 1 < N; m++)
        for (size_t r = 1; r <= m + 1; r++) {
            mpz_ptr d = x->stir[(m + 1) * N + r];
            mpz_set(d, x->stir[m * N + (r - 1)]);
            if (r <= m) mpz_submul_ui(d, x->stir[m * N + r], (unsigned long)m);
        }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Terms of a level                                                   *
 * ------------------------------------------------------------------ */

static bool nsb_block(NsbCtx* x, const NsbLevel* P, const mpz_t A, const mpz_t B,
                      mpfr_t out, double* err);

/* Start of block k of a structure (b, c): b^(k-1) + c.  False past the index
 * budget. */
static bool nsb_block_start(NsbCtx* x, const mpz_t b, long c, const mpz_t k, mpz_t out) {
    if (mpz_sgn(k) <= 0 || !mpz_fits_ulong_p(k)) { x->budget_hit = true; return false; }
    unsigned long e = mpz_get_ui(k) - 1;
    double bits = (double)e * nsb_log2_z(b);
    if (bits > (double)x->max_index_bits) { x->budget_hit = true; return false; }
    mpz_pow_ui(out, b, e);
    if (c >= 0) mpz_add_ui(out, out, (unsigned long)c);
    else        mpz_sub_ui(out, out, (unsigned long)(-c));
    return true;
}

/* Term k of level L at the working precision; *err accumulates its error. */
static bool nsb_term(NsbCtx* x, const NsbLevel* L, const mpz_t k, mpfr_t out, double* err) {
    if (nsb_over_budget(x)) return false;
    if (L->depth == 0) {
        if (mpz_sizeinbase(k, 2) > x->max_index_bits) { x->budget_hit = true; return false; }
        x->nevals++;
        if (!x->f(x->ud, k, out) || !mpfr_number_p(out)) { x->bad = true; return false; }
        return true;
    }
    mpz_t A, B, k1;
    mpz_inits(A, B, k1, (mpz_ptr)0);
    mpz_add_ui(k1, k, 1);
    /* A block sum is wanted to the accuracy its consumer needs: this level's
     * own, or less when it is computed as interpolation data for a coarser
     * level (x->gq, set by the caller). */
    long gsave = x->gq;
    if (L->tq < x->gq) x->gq = L->tq;
    bool ok = nsb_block_start(x, L->b, L->c, k, A)
           && nsb_block_start(x, L->b, L->c, k1, B)
           && nsb_block(x, L->parent, A, B, out, err);
    x->gq = gsave;
    if (ok && nsb_debug() > 1 && L->depth == 1)
        mpfr_fprintf(stderr, "[nsb]   B[%Zd] = %.30Rg evals %ld\n", k, out, x->nevals);
    mpz_clears(A, B, k1, (mpz_ptr)0);
    return ok;
}

/* Sum_{n=A}^{B-1} term(P, n) one term at a time. */
static bool nsb_direct(NsbCtx* x, const NsbLevel* P, const mpz_t A, const mpz_t B,
                       mpfr_t out, double* err) {
    mpfr_t t;
    mpfr_init2(t, x->p);
    mpz_t n;
    mpz_init_set(n, A);
    mpfr_set_ui(out, 0, MPFR_RNDN);
    bool ok = true;
    for (; mpz_cmp(n, B) < 0; mpz_add_ui(n, n, 1)) {
        if (!nsb_term(x, P, n, t, err)) { ok = false; break; }
        mpfr_add(out, out, t, MPFR_RNDN);
    }
    mpz_clear(n);
    mpfr_clear(t);
    return ok;
}

/* The smooth continuation of P's terms at real X in [A, E] (all integers of
 * [A, E] lie in one smooth block): a q-point Lagrange interpolant on the
 * integers nearest X, with q the smallest width whose error bound
 * prod|X - x_i| / X^q (the bound for a function smooth on the scale of its
 * argument, c/x and its kin) is below the target.  *rel receives that bound. */
static bool nsb_interp(NsbCtx* x, const NsbLevel* P, const mpz_t A, const mpz_t E,
                       const mpfr_t X, mpfr_t out, double* err) {
    double lx = nsb_log2(X);
    double goal = -(double)(x->gq + 12);
    mpz_t fl, s0, span, n;
    mpz_inits(fl, s0, span, n, (mpz_ptr)0);
    mpfr_get_z(fl, X, MPFR_RNDD);
    mpz_sub(span, E, A);                         /* E - A + 1 integers available */
    mpfr_t t;
    mpfr_init2(t, x->p);
    int q = 1;
    double bnd = 0.0;
    for (q = 1; q <= NSB_QMAX; q++) {
        if (mpz_cmp_ui(span, (unsigned long)(q - 1)) < 0) { q--; break; }
        mpz_sub_ui(s0, fl, (unsigned long)((q - 1) / 2));
        if (mpz_cmp(s0, A) < 0) mpz_set(s0, A);
        mpz_add_ui(n, s0, (unsigned long)(q - 1));
        if (mpz_cmp(n, E) > 0) { mpz_sub_ui(s0, E, (unsigned long)(q - 1)); }
        mpfr_sub_z(t, X, s0, MPFR_RNDN);
        double tt = mpfr_get_d(t, MPFR_RNDN);
        bnd = -(double)q * lx;
        for (int i = 0; i < q; i++) bnd += log2(fabs(tt - (double)i));
        if (bnd < goal) break;
    }
    if (q > NSB_QMAX) q = NSB_QMAX;
    if (q < 1) q = 1;
    /* Recompute the stencil for the chosen q (the loop may have run past it). */
    mpz_sub_ui(s0, fl, (unsigned long)((q - 1) / 2));
    if (mpz_cmp(s0, A) < 0) mpz_set(s0, A);
    mpz_add_ui(n, s0, (unsigned long)(q - 1));
    if (mpz_cmp(n, E) > 0) mpz_sub_ui(s0, E, (unsigned long)(q - 1));
    mpfr_sub_z(t, X, s0, MPFR_RNDN);

    bool ok = true;
    if (mpfr_integer_p(t) && mpfr_sgn(t) >= 0 && mpfr_cmp_ui(t, (unsigned long)q) < 0) {
        mpz_add_ui(n, s0, mpfr_get_ui(t, MPFR_RNDN));
        ok = nsb_term(x, P, n, out, err);
    } else {
        /* Barycentric form: P(t) = l(t) sum_i w_i y_i / (t - i),
         * w_i = (-1)^(q-1-i) / (i! (q-1-i)!), l(t) = prod (t - j). */
        mpfr_t y, acc, l, d, w, f1;
        mpfr_inits2(x->p, y, acc, l, d, w, f1, (mpfr_ptr)0);
        mpfr_set_ui(acc, 0, MPFR_RNDN);
        mpfr_set_ui(l, 1, MPFR_RNDN);
        for (int i = 0; i < q && ok; i++) {
            mpz_add_ui(n, s0, (unsigned long)i);
            if (!nsb_term(x, P, n, y, err)) { ok = false; break; }
            mpfr_sub_ui(d, t, (unsigned long)i, MPFR_RNDN);
            mpfr_mul(l, l, d, MPFR_RNDN);
            mpfr_fac_ui(w, (unsigned long)i, MPFR_RNDN);
            mpfr_fac_ui(f1, (unsigned long)(q - 1 - i), MPFR_RNDN);
            mpfr_mul(w, w, f1, MPFR_RNDN);
            mpfr_mul(w, w, d, MPFR_RNDN);
            mpfr_div(y, y, w, MPFR_RNDN);
            if ((q - 1 - i) & 1) mpfr_sub(acc, acc, y, MPFR_RNDN);
            else                 mpfr_add(acc, acc, y, MPFR_RNDN);
        }
        if (ok) mpfr_mul(out, acc, l, MPFR_RNDN);
        mpfr_clears(y, acc, l, d, w, f1, (mpfr_ptr)0);
    }
    if (ok && bnd > goal) *err += nsb_absd(out) * exp2(bnd);
    mpfr_clear(t);
    mpz_clears(fl, s0, span, n, (mpz_ptr)0);
    return ok;
}

/* Derivative of order r at the stencil's first point from the q samples
 * y[0..q-1] at unit spacing: D^r = r! sum_{m>=r} s(m,r) Delta^m / m!
 * (forward, sgn = +1), or with |s(m,r)| and backward differences when the
 * samples run away from the point (sgn = -1, y[i] = G(E - i)). */
static void nsb_deriv(NsbCtx* x, mpfr_t* y, int q, int r, int sgn, mpfr_t out) {
    mpfr_t dm, term, cz;
    mpfr_inits2(x->p, dm, term, cz, (mpfr_ptr)0);
    mpz_t bin;
    mpz_init(bin);
    mpfr_set_ui(out, 0, MPFR_RNDN);
    for (int m = r; m < q; m++) {
        /* forward: Delta^m y0 = sum_i (-1)^(m-i) C(m,i) y_i
         * backward: nabla^m yE = sum_i (-1)^i C(m,i) y_i (y_i = G(E - i)) */
        mpfr_set_ui(dm, 0, MPFR_RNDN);
        for (int i = 0; i <= m; i++) {
            mpz_bin_uiui(bin, (unsigned long)m, (unsigned long)i);
            mpfr_mul_z(term, y[i], bin, MPFR_RNDN);
            bool neg = (sgn > 0) ? ((m - i) & 1) : (i & 1);
            if (neg) mpfr_sub(dm, dm, term, MPFR_RNDN);
            else     mpfr_add(dm, dm, term, MPFR_RNDN);
        }
        mpz_set(bin, x->stir[(size_t)m * NSB_QE_MAX + (size_t)r]);
        if (sgn < 0) mpz_abs(bin, bin);
        mpfr_mul_z(term, dm, bin, MPFR_RNDN);
        mpfr_fac_ui(cz, (unsigned long)m, MPFR_RNDN);
        mpfr_div(term, term, cz, MPFR_RNDN);
        mpfr_add(out, out, term, MPFR_RNDN);
    }
    mpfr_fac_ui(cz, (unsigned long)r, MPFR_RNDN);
    mpfr_mul(out, out, cz, MPFR_RNDN);
    mpz_clear(bin);
    mpfr_clears(dm, term, cz, (mpfr_ptr)0);
}

/* Euler–Maclaurin for Sum_{n=A}^{E} G(n), G smooth on [A, E], A >= T:
 *   int_A^E G + (G(A) + G(E))/2 + sum_j B_2j/(2j)! (G^(2j-1)(E) - G^(2j-1)(A)),
 * the integral in u = log n. */
static bool nsb_em(NsbCtx* x, const NsbLevel* P, const mpz_t A, const mpz_t E,
                   mpfr_t out, double* err) {
    double lA = nsb_log2_z(A);
    double goal = -(double)(x->gq + 12);
    /* Correction j relative to the block, ~ |B_2j| / (2j A^(2j+1)). */
    int J = 0;
    for (int j = 1; j <= NSB_JMAX; j++) {
        double lb = 1.0 + lgamma(2.0 * j + 1.0) / log(2.0) - 2.0 * j * log2(2.0 * M_PI)
                    - log2(2.0 * j) - (2.0 * j + 1.0) * lA;
        if (lb < goal) break;
        J = j;
    }
    int qe = J > 0 ? 2 * J + 6 : 1;
    if (qe > NSB_QE_MAX) qe = NSB_QE_MAX;

    bool ok = true;
    mpfr_t* yA = malloc(sizeof(mpfr_t) * (size_t)qe);
    mpfr_t* yE = malloc(sizeof(mpfr_t) * (size_t)qe);
    if (!yA || !yE) { free(yA); free(yE); return false; }
    for (int i = 0; i < qe; i++) { mpfr_init2(yA[i], x->p); mpfr_init2(yE[i], x->p); }
    mpz_t n;
    mpz_init(n);
    for (int i = 0; i < qe && ok; i++) {
        mpz_add_ui(n, A, (unsigned long)i);
        ok = nsb_term(x, P, n, yA[i], err);
        if (ok) { mpz_sub_ui(n, E, (unsigned long)i); ok = nsb_term(x, P, n, yE[i], err); }
    }

    mpfr_t acc, dA, dE, tmp, uA, uE, mid, half, X, g, Iq, Iprev;
    mpfr_inits2(x->p, acc, dA, dE, tmp, uA, uE, mid, half, X, g, Iq, Iprev, (mpfr_ptr)0);
    if (ok) {
        /* (G(A) + G(E)) / 2 and the derivative corrections. */
        mpfr_add(acc, yA[0], yE[0], MPFR_RNDN);
        mpfr_div_ui(acc, acc, 2, MPFR_RNDN);
        for (int j = 1; j <= J; j++) {
            nsb_deriv(x, yA, qe, 2 * j - 1, +1, dA);
            nsb_deriv(x, yE, qe, 2 * j - 1, -1, dE);
            mpfr_sub(tmp, dE, dA, MPFR_RNDN);
            mpfr_mul(tmp, tmp, x->bern[j], MPFR_RNDN);
            mpfr_add(acc, acc, tmp, MPFR_RNDN);
            /* Truncation: the next term is the last one times about
             * (2J+1)(2J+2) / (2 pi A)^2 (B_2j/(2j)! ~ (2 pi)^-2j, G^(r) ~ r! G/A^r). */
            if (j == J) {
                double r = (2.0 * J + 1.0) * (2.0 * J + 2.0)
                         / exp2(2.0 * (lA + log2(2.0 * M_PI)));
                *err += nsb_absd(tmp) * 4.0 * r;
            }
        }
        /* The integral in u = log n by Gauss–Legendre rules of doubling size. */
        mpfr_set_z(tmp, A, MPFR_RNDN); mpfr_log(uA, tmp, MPFR_RNDN);
        mpfr_set_z(tmp, E, MPFR_RNDN); mpfr_log(uE, tmp, MPFR_RNDN);
        mpfr_add(mid, uA, uE, MPFR_RNDN); mpfr_div_ui(mid, mid, 2, MPFR_RNDN);
        mpfr_sub(half, uE, uA, MPFR_RNDN); mpfr_div_ui(half, half, 2, MPFR_RNDN);
        bool have = false;
        double diff = INFINITY, ierr = 0.0;
        int nrule = 0;
        for (int r = 0; r < NSB_GL_RULES && ok; r++) {
            nrule = r;
            if (!nsb_gl_build(x, r)) { ok = false; break; }
            int nn = x->gl_n[r];
            ierr = 0.0;
            mpfr_set_ui(Iq, 0, MPFR_RNDN);
            for (int i = 0; i < nn && ok; i++) {
                mpfr_mul(tmp, half, x->gl_x[r][i], MPFR_RNDN);
                mpfr_add(tmp, mid, tmp, MPFR_RNDN);
                mpfr_exp(X, tmp, MPFR_RNDN);
                /* Rounding of exp may step outside [A, E] at the ends. */
                if (mpfr_cmp_z(X, A) < 0) mpfr_set_z(X, A, MPFR_RNDN);
                if (mpfr_cmp_z(X, E) > 0) mpfr_set_z(X, E, MPFR_RNDN);
                double enode = 0.0;          /* absolute error of G(X) */
                ok = nsb_interp(x, P, A, E, X, g, &enode);
                if (!ok) break;
                mpfr_mul(g, g, X, MPFR_RNDN);
                mpfr_mul(g, g, x->gl_w[r][i], MPFR_RNDN);
                mpfr_add(Iq, Iq, g, MPFR_RNDN);
                ierr += enode * mpfr_get_d(X, MPFR_RNDU) * mpfr_get_d(x->gl_w[r][i], MPFR_RNDU);
            }
            if (!ok) break;
            mpfr_mul(Iq, Iq, half, MPFR_RNDN);
            ierr *= nsb_absd(half);
            if (have) {
                mpfr_sub(tmp, Iq, Iprev, MPFR_RNDN);
                diff = nsb_absd(tmp);
                /* Converged once the finer rule's QUADPACK estimate (below)
                 * meets the goal: |I| (200 d/|I|)^1.5 <= 2^goal |S|.  (Asking d
                 * itself for that would chase the interpolants' noise floor
                 * through every rule.) */
                double scale = nsb_absd(Iq) + nsb_absd(acc);
                if (diff <= exp2(goal / 1.5) * scale / 200.0) break;
            }
            mpfr_set(Iprev, Iq, MPFR_RNDN);
            have = true;
        }
        if (ok) {
            /* The doubled rule's difference d bounds the coarser rule's error;
             * for an analytic integrand the finer rule's is ~ d^2/|I|.  The
             * QUADPACK form |I| min(1, (200 d/|I|)^1.5) sits between the two,
             * with a floor at the target resolution. */
            double aI = nsb_absd(Iq);
            double qe = isfinite(diff) ? (aI > 0.0 ? aI * fmin(1.0, pow(200.0 * diff / aI, 1.5))
                                                    : diff)
                                       : aI;
            qe = fmax(qe, aI * exp2(goal));
            if (nsb_debug() > 2 || (nsb_debug() > 1 && P->depth > 0))
                gmp_fprintf(stderr, "[nsb]     em depth %d A=%Zd J=%d rule %d diff %g qe %g ierr %g |I| %g\n",
                            P->depth, A, J, nrule, diff, qe, ierr, aI);
            *err += qe + (isfinite(ierr) ? ierr : 0.0);
            mpfr_add(out, Iq, acc, MPFR_RNDN);
        }
    }
    mpfr_clears(acc, dA, dE, tmp, uA, uE, mid, half, X, g, Iq, Iprev, (mpfr_ptr)0);
    for (int i = 0; i < qe; i++) { mpfr_clear(yA[i]); mpfr_clear(yE[i]); }
    free(yA); free(yE);
    mpz_clear(n);
    return ok;
}

/* The narrowest centred stencil that interpolates at x to the goal, by the
 * bound ((q/2)!)^2 / x^q of nsb_interp; NSB_QMAX + 1 when none does. */
static int nsb_q_needed(const NsbCtx* x, double lx) {
    double goal = -(double)(x->gq + 12);
    for (int q = 1; q <= NSB_QMAX; q++)
        if (2.0 * lgamma(q / 2.0 + 1.0) / log(2.0) - (double)q * lx < goal) return q;
    return NSB_QMAX + 1;
}

/* Sum_{n=A}^{B-1} of P's terms, which are smooth on [A, B).  Term by term
 * below NSB_EM_FLOOR, below the point where a stencil of at most NSB_QMAX
 * reaches the goal, and wherever that is cheaper than Euler–Maclaurin's ~40
 * interpolated samples per stencil point; Euler–Maclaurin beyond. */
#define NSB_EM_FLOOR 64
static bool nsb_block(NsbCtx* x, const NsbLevel* P, const mpz_t A, const mpz_t B,
                      mpfr_t out, double* err) {
    if (mpz_cmp(A, B) >= 0) { mpfr_set_ui(out, 0, MPFR_RNDN); return true; }
    mpz_t A2, E, len;
    mpz_inits(A2, E, len, (mpz_ptr)0);
    mpz_set(A2, A);
    if (mpz_cmp_ui(A2, NSB_EM_FLOOR) < 0) mpz_set_ui(A2, NSB_EM_FLOOR);
    int q = nsb_q_needed(x, nsb_log2_z(A2));
    if (q > NSB_QMAX) {
        /* first x where the widest stencil suffices */
        double lx = ((double)(x->gq + 12) + 2.0 * lgamma(NSB_QMAX / 2.0 + 1.0) / log(2.0))
                    / (double)NSB_QMAX + 0.01;
        mpz_t t;
        mpz_init(t);
        mpz_set_d(t, ceil(exp2(lx)));
        if (mpz_cmp(t, A2) > 0) mpz_set(A2, t);
        mpz_clear(t);
        q = nsb_q_needed(x, nsb_log2_z(A2));
    }
    if (mpz_cmp(A2, B) > 0) mpz_set(A2, B);
    mpz_sub(len, B, A2);
    double em_cost = 40.0 * (double)q + 2.0 * NSB_QE_MAX;
    if (mpz_cmp_d(len, em_cost) <= 0 || mpz_cmp_ui(len, NSB_DIRECT_MAX) <= 0) mpz_set(A2, B);
    bool ok = nsb_direct(x, P, A, A2, out, err);
    if (ok && mpz_cmp(A2, B) < 0) {
        mpfr_t s;
        mpfr_init2(s, x->p);
        mpz_sub_ui(E, B, 1);
        ok = nsb_em(x, P, A2, E, s, err);
        if (ok) mpfr_add(out, out, s, MPFR_RNDN);
        mpfr_clear(s);
    }
    mpz_clears(A2, E, len, (mpz_ptr)0);
    return ok;
}

/* ------------------------------------------------------------------ *
 *  Structure of a level's head                                        *
 * ------------------------------------------------------------------ */

typedef enum { NSB_SMOOTH, NSB_BLOCKED, NSB_MORE, NSB_ROUGH } NsbShape;

/* Is u = b^j for some j >= 1? */
static bool nsb_is_power(unsigned long u, unsigned long b) {
    if (b < 2 || u < b) return false;
    while (u % b == 0) u /= b;
    return u == 1;
}

/* Classify the terms a[0..n-1] at indices k0.. (k0 small, fits a long).
 * A breakpoint at m: log|a| follows a smooth law up to m-1 and another from
 * m on, each side predicting its own end (cubic extrapolation, error
 * < NSB_FLAG_SMALL) but not the other's (error > NSB_FLAG_BIG) — or the sign
 * flips between two runs of constant sign.  The breakpoints must be exactly
 * the points b^j + c (integer b >= 2, |c| <= 2) the head can resolve, at
 * least two of them; any irregularity off those points is ROUGH. */
static NsbShape nsb_classify(mpfr_t* a, int n, long k0, unsigned long* bout, long* cout) {
    if (n < 16) return NSB_MORE;
    double* l = malloc(sizeof(double) * (size_t)n);
    int* sg = malloc(sizeof(int) * (size_t)n);
    long* F = malloc(sizeof(long) * (size_t)n);
    if (!l || !sg || !F) { free(l); free(sg); free(F); return NSB_ROUGH; }
    double lmax = 0.0;
    NsbShape shape = NSB_ROUGH;
    int nf = 0;
    for (int i = 0; i < n; i++) {
        if (mpfr_zero_p(a[i])) goto done;
        l[i] = nsb_log2(a[i]) * log(2.0);
        sg[i] = mpfr_sgn(a[i]);
        if (fabs(l[i]) > lmax) lmax = fabs(l[i]);
    }
    {
        double floor_e = 1e-11 * (1.0 + lmax) * 16.0;
#define EL(i) fabs(l[i] - (4*l[(i)-1] - 6*l[(i)-2] + 4*l[(i)-3] - l[(i)-4]))
#define ER(i) fabs(l[i] - (4*l[(i)+1] - 6*l[(i)+2] + 4*l[(i)+3] - l[(i)+4]))
        bool* isflag = calloc((size_t)n, sizeof(bool));
        if (!isflag) goto done;
        for (int i = 5; i + 5 < n; i++) {
            bool sides = ER(i) < NSB_FLAG_SMALL && EL(i - 1) < NSB_FLAG_SMALL;
            bool mag = EL(i) > NSB_FLAG_BIG && ER(i - 1) > NSB_FLAG_BIG;
            bool sgn = sg[i] != sg[i - 1];
            if (sgn) {
                for (int j = i - 4; j < i - 1; j++) if (sg[j] != sg[j + 1]) sgn = false;
                for (int j = i; j < i + 3; j++) if (sg[j] != sg[j + 1]) sgn = false;
            }
            if (sides && (mag || sgn)) { isflag[i] = true; F[nf++] = k0 + i; }
        }
        /* Roughness off the breakpoints: a smooth law has fourth differences of
         * log|a| ~ 1/m^4; anything far above that (prime gaps, divisor counts)
         * is not a law Levin or a block quadrature can use. */
        for (int i = 20; i + 5 < n; i++) {
            bool near = false;
            for (int j = i - 5; j <= i + 5; j++) if (j >= 0 && j < n && isflag[j]) near = true;
            if (near) continue;
            double m = (double)(k0 + i);
            double e = fmin(EL(i), ER(i));
            if (e > floor_e && e > 1e3 / (m * m * m * m) + floor_e) { free(isflag); goto done; }
        }
        free(isflag);
#undef EL
#undef ER
    }
    /* Signs: Levin's model wants a constant or a strictly alternating sign,
     * and a block sum by Euler–Maclaurin one sign per block (checked against
     * the breakpoints below).  (-1)^n/f(BitLength[n]) is neither: ROUGH. */
    int nchg = 0;
    for (int i = 1; i < n; i++) if (sg[i] != sg[i - 1]) nchg++;
    if (nf == 0) {
        if (nchg == 0 || nchg == n - 1) shape = NSB_SMOOTH;
        goto done;
    }
    if (nchg >= n / 2) goto done;
    if (nf < 2) { shape = NSB_MORE; goto done; }
    {
        static const long cs[5] = { 0, -1, 1, -2, 2 };
        long m1 = F[nf - 2], m2 = F[nf - 1];
        for (int ci = 0; ci < 5; ci++) {
            long c = cs[ci];
            long u1 = m1 - c, u2 = m2 - c;
            if (u1 < 2 || u2 % u1 != 0) continue;
            unsigned long b = (unsigned long)(u2 / u1);
            if (!nsb_is_power((unsigned long)u1, b) || !nsb_is_power((unsigned long)u2, b)) continue;
            /* every flag must be a b^j + c, and so every sign change ... */
            bool good = true;
            for (int i = 0; i < nf && good; i++)
                if (F[i] - c < 2 || !nsb_is_power((unsigned long)(F[i] - c), b)) good = false;
            for (int i = 1; i < n && good; i++) {
                long m = k0 + i;
                if (sg[i] != sg[i - 1] && (m - c < 2 || !nsb_is_power((unsigned long)(m - c), b)))
                    good = false;
            }
            /* ... and every resolvable b^j + c a flag (both neighbouring blocks
             * long enough for the 5-point stencils, inside the window). */
            int need = 0;
            for (unsigned long pw = b; good && (long)pw + c <= k0 + n - 6; ) {
                long m = (long)pw + c;
                if (m >= k0 + 5 && pw - pw / b >= 5) {
                    need++;
                    bool found = false;
                    for (int i = 0; i < nf; i++) if (F[i] == m) found = true;
                    if (!found) good = false;
                }
                if (pw > (unsigned long)(LONG_MAX / 2) / b) break;
                pw *= b;
            }
            if (good && need >= 2) { *bout = b; *cout = c; shape = NSB_BLOCKED; goto done; }
        }
        shape = NSB_MORE;   /* breakpoints, not (yet) a verified geometric pattern */
    }
done:
    free(l); free(sg); free(F);
    return shape;
}

/* ------------------------------------------------------------------ *
 *  One level                                                          *
 * ------------------------------------------------------------------ */

/* Levin length for the target. */
static int nsb_levin_len(const NsbCtx* x) {
    return (int)((double)x->tb * 0.30103 * 1.3) + 10;
}

/* Levin's u-transform on the partial sums of a[s..s+L-1] (index k0+s on),
 * plus the direct head a[0..s-1].  *e gets the transform's own last step and
 * *gam its stability factor sum_j |w_j/omega_j| / |sum_j w_j/omega_j| — the
 * factor by which an error in the partial sums reaches the result, which for
 * a logarithmically convergent series grows steeply with the start index
 * (high-order differences of a law smooth on the scale k0 + s). */
static bool nsb_levin_at(NsbCtx* x, mpfr_t* a, int s, int L, long k0, mpfr_t out,
                         double* e, double* gam) {
    mpfr_t* S = malloc(sizeof(mpfr_t) * (size_t)L);
    mpfr_t* Z = malloc(sizeof(mpfr_t) * (size_t)L);
    if (!S || !Z) { free(S); free(Z); return false; }
    for (int i = 0; i < L; i++) {
        mpfr_init2(S[i], x->p); mpfr_init2(Z[i], x->p);
        mpfr_set_ui(Z[i], 0, MPFR_RNDN);
        if (i == 0) mpfr_set(S[0], a[s], MPFR_RNDN);
        else        mpfr_add(S[i], S[i - 1], a[s + i], MPFR_RNDN);
    }
    mpfr_t ri, ii;
    mpfr_inits2(x->p, ri, ii, (mpfr_ptr)0);
    double st = INFINITY;
    bool fin = false;
    double beta = (double)(k0 + s);
    bool ok = seqaccel_levin_mpfr(S, Z, L, SEQACCEL_LEVIN_U, beta, (long)x->p,
                                  ri, ii, &st, &fin) && fin;
    if (ok) {
        for (int i = 0; i < s; i++) mpfr_add(ri, ri, a[i], MPFR_RNDN);
        mpfr_set(out, ri, MPFR_RNDN);
        *e = st;
        /* Stability factor of the degree k = L-2 transform from base n = 1
         * (seqaccel's convention): w_j = (-1)^j C(k,j) ((beta+1+j)/(beta+1+k))^(k-1),
         * omega_{1+j} = (beta+1+j) a_{1+j}. */
        int k = L - 2;
        mpfr_t w, q, sa, sd;
        mpfr_inits2(x->p, w, q, sa, sd, (mpfr_ptr)0);
        mpz_t bz;
        mpz_init(bz);
        mpfr_set_ui(sa, 0, MPFR_RNDN);
        mpfr_set_ui(sd, 0, MPFR_RNDN);
        for (int j = 0; j <= k; j++) {
            mpfr_set_d(q, (beta + 1.0 + j) / (beta + 1.0 + k), MPFR_RNDN);
            mpfr_pow_ui(w, q, (unsigned long)(k - 1), MPFR_RNDN);
            mpz_bin_uiui(bz, (unsigned long)k, (unsigned long)j);
            mpfr_mul_z(w, w, bz, MPFR_RNDN);
            mpfr_mul_d(q, a[s + 1 + j], beta + 1.0 + j, MPFR_RNDN);
            mpfr_div(w, w, q, MPFR_RNDN);
            if (j & 1) mpfr_sub(sd, sd, w, MPFR_RNDN); else mpfr_add(sd, sd, w, MPFR_RNDN);
            mpfr_abs(w, w, MPFR_RNDN);
            mpfr_add(sa, sa, w, MPFR_RNDN);
        }
        *gam = mpfr_zero_p(sd) ? INFINITY : exp2(nsb_log2(sa) - nsb_log2(sd));
        mpfr_clears(w, q, sa, sd, (mpfr_ptr)0);
        mpz_clear(bz);
    }
    mpfr_clears(ri, ii, (mpfr_ptr)0);
    for (int i = 0; i < L; i++) { mpfr_clear(S[i]); mpfr_clear(Z[i]); }
    free(S); free(Z);
    return ok;
}

/* Levin on the head a[0..n-1], from the best of several starts s = 0, 2, 4,
 * 8, ...: the first terms of a series often carry a transient its asymptotic
 * law does not (a block sum's psi corrections ~ b^-k), which Levin's model
 * cannot absorb, so the head below s is summed directly — while later starts
 * amplify the data error more.  Each start's error is the larger of the
 * transform's last step and its stability factor times `dS`, the error of the
 * partial sums; the start with the smallest is taken, and then widened to its
 * disagreement with any other start the two error bars cannot reconcile
 * (within a factor of ten). */
#define NSB_LEVIN_STARTS 9
static bool nsb_levin(NsbCtx* x, mpfr_t* a, int n, long k0, double dS, int s0,
                      mpfr_t out, double* err, int* want) {
    int L = nsb_levin_len(x);
    if (L > n) L = n;
    *want = 0;
    if (L < 8) return false;
    /* s0 > 0: the first start past the level's known transient (see
     * nsb_transient_start); the doubling ladder from there. */
    int st[NSB_LEVIN_STARTS] = { 0, 2, 4, 8, 16, 32, 64, 128, 256 };
    if (s0 > 0) {
        st[0] = 0;
        for (int i = 1; i < NSB_LEVIN_STARTS; i++) st[i] = s0 + (i == 1 ? 0 : 4 << (i - 2));
    }
    mpfr_t E[NSB_LEVIN_STARTS];
    double es[NSB_LEVIN_STARTS];
    int m = 0;
    /* At each start, four degrees: a late start (large beta) is better served
     * by a lower one — the transform's high-order differences of a law smooth
     * on the scale beta lose to its conditioning (Sum 1/k^2 from k = 248:
     * 1e-72 with 40 terms, 1e-26 with 88). */
    mpfr_t Ec;
    mpfr_init2(Ec, x->p);
    for (int i = 0; i < NSB_LEVIN_STARTS && st[i] + L <= n; i++) {
        static const int num[4] = { 3, 2, 3, 1 }, den[4] = { 3, 3, 6, 3 };
        bool any = false;
        mpfr_init2(E[m], x->p);
        es[m] = INFINITY;
        for (int v = 0; v < 4; v++) {
            int Lc = L * num[v] / den[v];
            if (Lc < 10) continue;
            double e = INFINITY, g = INFINITY;
            if (!nsb_levin_at(x, a, st[i], Lc, k0, Ec, &e, &g)) continue;
            /* The last step under-reports by a small factor as a rule. */
            double ec = fmax(4.0 * e, g * dS);
            if (!isfinite(ec)) ec = INFINITY;
            if (nsb_debug() > 1)
                mpfr_fprintf(stderr, "[nsb]   levin start %d len %d: %.25Rg step %g gamma %g err %g\n",
                             st[i], Lc, Ec, e, g, ec);
            if (!any || ec < es[m]) { mpfr_set(E[m], Ec, MPFR_RNDN); es[m] = ec; any = true; }
        }
        if (!any) { mpfr_clear(E[m]); break; }
        m++;
    }
    mpfr_clear(Ec);
    if (m == 0) return false;
    int best = 0;
    for (int i = 1; i < m; i++) if (es[i] < es[best]) best = i;
    double beste = es[best];
    if (m == 1) beste = fmax(beste, nsb_absd(a[L - 1]) * (double)L);
    mpfr_t d;
    mpfr_init2(d, x->p);
    for (int i = 0; i < m; i++) {
        if (i == best || !isfinite(es[i])) continue;
        mpfr_sub(d, E[i], E[best], MPFR_RNDN);
        double dd = nsb_absd(d);
        /* A transform's last step under-reports its error by a few times as a
         * matter of course; ten times the combined bars is a real conflict. */
        if (dd > 10.0 * (es[i] + es[best]) && dd > beste) beste = dd;
    }
    /* The latest start has nothing after it to cross-check it: ask for the
     * terms the next start would need, so the caller can fetch them. */
    if (best == m - 1 && m < NSB_LEVIN_STARTS) *want = st[m] + L;
    mpfr_set(out, E[best], MPFR_RNDN);
    *err += beste;
    mpfr_clear(d);
    for (int i = 0; i < m; i++) mpfr_clear(E[i]);
    return true;
}

/* The terms of a level have stopped mattering: the last few each fall below
 * the target relative to the sum, decreasing at least geometrically. */
static bool nsb_settled(NsbCtx* x, mpfr_t* a, int n, const mpfr_t S, double* tail) {
    if (n < 6) return false;
    double lS = nsb_log2(S);
    if (!isfinite(lS)) return false;
    double lr = -INFINITY;
    for (int i = n - 4; i < n; i++) {
        double li = nsb_log2(a[i]), lp = nsb_log2(a[i - 1]);
        if (li > lS - (double)(x->tb + 10)) return false;
        if (li - lp > log2(0.9)) return false;
        if (li - lp > lr) lr = li - lp;
    }
    double r = exp2(lr);
    *tail = exp2(nsb_log2(a[n - 1])) * r / (1.0 - r);
    return true;
}

/* A bound on the unsummed remainder Sum_{k > k0+n-1} a_k of a head that
 * could not be summed to convergence, from its last terms: |a| when the last
 * two alternate in sign and shrink; otherwise the slowest power law k^-p
 * through consecutive pairs of the last five, whose tail is |a| k/(p-1),
 * charged twice — and INFINITY when p <= 1.2: a
 * remainder decaying that slowly (or not at all, as far as the head shows)
 * has no bound the head can give. */
static double nsb_powtail(mpfr_t* a, int n, long k0) {
    if (n < 2) return INFINITY;
    double l1 = nsb_log2(a[n - 1]), l0 = nsb_log2(a[n - 2]);
    if (!isfinite(l1)) return 0.0;
    if (mpfr_sgn(a[n - 1]) != mpfr_sgn(a[n - 2]))
        return l1 < l0 ? exp2(l1) : INFINITY;
    /* The slowest local exponent over the last few pairs: one pair can
     * straddle a breakpoint of a structure the head is too short to see. */
    double p = INFINITY;
    for (int i = n - 1; i >= 1 && i >= n - 4; i--) {
        if (mpfr_sgn(a[i]) != mpfr_sgn(a[i - 1])) return INFINITY;
        double k = (double)(k0 + i);
        double pi = (nsb_log2(a[i - 1]) - nsb_log2(a[i])) / log2(k / (k - 1.0));
        if (pi < p) p = pi;
    }
    if (!(p > 1.2)) return INFINITY;
    return 2.0 * exp2(l1) * (double)(k0 + n - 1) / (p - 1.0);
}

/* A level of block sums (depth >= 1) carries, besides the law it converges
 * by, the Euler–Maclaurin endpoint terms of each block, which fall off like
 * b^-k: real terms of the series, but no part of the algebraic model Levin
 * extrapolates, and amplified by it.  The first index past which they are
 * below the target, relative: (tb + 16) / log2 b; as a start offset from k0. */
static int nsb_transient_start(const NsbCtx* x, const NsbLevel* L, long k0) {
    if (L->depth == 0) return 0;
    double lb = nsb_log2_z(L->b);
    if (!(lb > 0.0)) return 0;
    long k = (long)ceil(((double)x->tb + 16.0) / lb);
    return k > k0 ? (int)(k - k0) : 0;
}

/* Extend the head a[0..*n-1] of level L (indices k0..) to `upto` terms,
 * accumulating S.  Stops early on a failed term (budget, or x->bad), or — when
 * `settle` — once the terms have stopped mattering (*settled, *tail). */
static void nsb_fetch(NsbCtx* x, const NsbLevel* L, mpfr_t* a, int* n, int upto, long k0,
                      mpfr_t S, double* lerr, bool settle, bool* settled, double* tail) {
    mpz_t k;
    mpz_init(k);
    while (*n < upto) {
        mpfr_init2(a[*n], x->p);
        mpz_set_si(k, k0 + *n);
        if (!nsb_term(x, L, k, a[*n], lerr)) { mpfr_clear(a[*n]); break; }
        mpfr_add(S, S, a[*n], MPFR_RNDN);
        (*n)++;
        if (settle && nsb_settled(x, a, *n, S, tail)) { *settled = true; break; }
    }
    mpz_clear(k);
}

/* Sum_{k >= k0} term(L, k). */
static bool nsb_sum_level(NsbCtx* x, const NsbLevel* L, long k0, mpfr_t out, double* err) {
    if (L->depth > x->max_level) x->max_level = L->depth;
    long gsave = x->gq;
    x->gq = L->tq;
    int s0 = nsb_transient_start(x, L, k0);
    /* Level 0 starts on a short head and widens to NSB_H0 when it shows
     * breakpoints, or when its terms are cheap enough that looking three
     * decades out for a late breakpoint costs nothing (an expensive smooth
     * summand, BesselJZero[0, n], need not pay for 1100 terms). */
    int H = (L->depth == 0) ? NSB_H0_FIRST : s0 + nsb_levin_len(x) + 66;
    int Hmax = (L->depth == 0) ? NSB_H0 : (L->depth == 1 ? 1200 : 160);
    if (H > Hmax) Hmax = H;
    int cap = Hmax;
    mpfr_t* a = malloc(sizeof(mpfr_t) * (size_t)cap);
    if (!a) return false;
    int n = 0;
    mpfr_t S;
    mpfr_init2(S, x->p);
    mpfr_set_ui(S, 0, MPFR_RNDN);
    bool ok = false, settled = false;
    double tail = 0.0, lerr = 0.0;
    NsbShape shape = NSB_MORE;
    unsigned long b = 0; long c = 0;
    for (;;) {
        int nbefore = n;
        double t0 = mth_mono_seconds();
        nsb_fetch(x, L, a, &n, H, k0, S, &lerr, L->depth > 0, &settled, &tail);
        double per = (n > nbefore) ? (mth_mono_seconds() - t0) / (double)(n - nbefore) : INFINITY;
        if (x->bad) goto out;
        if (settled) {
            mpfr_set(out, S, MPFR_RNDN);
            *err += lerr + tail + nsb_absd(S) * exp2(-(double)(x->tb + 8));
            ok = true;
            goto out;
        }
        shape = nsb_classify(a, n, k0, &b, &c);
        if (nsb_debug())
            fprintf(stderr, "[nsb] depth %d k0 %ld n %d shape %d b %lu c %ld budget %d\n",
                    L->depth, k0, n, (int)shape, b, c, (int)x->budget_hit);
        if (L->depth == 0) {
            bool widen = n == H && H < Hmax && !x->budget_hit
                      && (shape == NSB_MORE || (shape == NSB_SMOOTH && per < NSB_CHEAP_TERM_S));
            if (widen) { H = Hmax; continue; }
            if (shape == NSB_BLOCKED || shape == NSB_SMOOTH) break;
            goto out;                                      /* decline: ROUGH / unverified */
        }
        if (shape == NSB_BLOCKED || shape == NSB_SMOOTH) break;
        if (n < H || H >= Hmax || x->budget_hit) break;    /* cannot look further */
        H = 3 * H + 10;
        if (H > Hmax) H = Hmax;
    }
    if (shape == NSB_BLOCKED && L->depth + 1 <= NSB_MAX_DEPTH) {
        /* The head up to the first block start >= NSB_EM_FLOOR, then the
         * block series (nsb_block decides, per block, what is summed term by
         * term). */
        mpz_t bz, pw, m;
        mpz_inits(bz, pw, m, (mpz_ptr)0);
        mpz_set_ui(bz, b);
        long j = 1;
        mpz_set(pw, bz);
        for (;;) {                                     /* m = b^j + c, first >= T, > k0 */
            mpz_set(m, pw);
            if (c >= 0) mpz_add_ui(m, m, (unsigned long)c); else mpz_sub_ui(m, m, (unsigned long)(-c));
            if (mpz_cmp_ui(m, NSB_EM_FLOOR) >= 0 && mpz_cmp_si(m, k0) > 0) break;
            mpz_mul(pw, pw, bz);
            j++;
        }
        /* head = window terms below m, then each block between the window end
         * and m (they are whole blocks: the window is contiguous from k0). */
        mpfr_t head, s;
        mpfr_inits2(x->p, head, s, (mpfr_ptr)0);
        mpfr_set_ui(head, 0, MPFR_RNDN);
        int i = 0;
        for (; i < n && mpz_cmp_si(m, k0 + i) > 0; i++) mpfr_add(head, head, a[i], MPFR_RNDN);
        bool hok = true;
        if (i == n && mpz_cmp_si(m, k0 + n) > 0) {
            /* blocks [from, next breakpoint) up to m */
            mpz_t from, nb, q;
            mpz_inits(from, nb, q, (mpz_ptr)0);
            mpz_set_si(from, k0 + n);
            while (hok && mpz_cmp(from, m) < 0) {
                /* next breakpoint > from */
                mpz_set(q, bz);
                for (;;) {
                    mpz_set(nb, q);
                    if (c >= 0) mpz_add_ui(nb, nb, (unsigned long)c); else mpz_sub_ui(nb, nb, (unsigned long)(-c));
                    if (mpz_cmp(nb, from) > 0) break;
                    mpz_mul(q, q, bz);
                }
                if (mpz_cmp(nb, m) > 0) mpz_set(nb, m);
                hok = nsb_block(x, L, from, nb, s, &lerr);
                if (hok) mpfr_add(head, head, s, MPFR_RNDN);
                mpz_set(from, nb);
            }
            mpz_clears(from, nb, q, (mpz_ptr)0);
        }
        if (hok) {
            NsbLevel child;
            child.depth = L->depth + 1;
            child.parent = L;
            mpz_init_set(child.b, bz);
            child.c = c;
            /* How far the index budget lets the block series run: its index i
             * needs this level's index b^i + c.  A series that cannot reach the
             * Levin length is closed by the short-head rule, which neither
             * amplifies its terms' errors nor, from a handful of terms, bounds
             * its remainder much below half the target: its terms are computed
             * to half the target bits (their quadrature error still charged). */
            {
                double r = L->lmax / nsb_log2_z(bz);
                child.lmax = r > 1.0 ? log2(r) : 0.0;
                double m = exp2(child.lmax) - (double)(j + 1) + 1.0;
                long half = x->tb / 2 > 53 ? x->tb / 2 : 53;
                child.tq = (m < (double)nsb_levin_len(x) + 8.0) ? half + 8 : x->tq;
            }
            double cerr = 0.0;
            bool cok = nsb_sum_level(x, &child, j + 1, s, &cerr);
            mpz_clear(child.b);
            if (nsb_debug())
                mpfr_fprintf(stderr, "[nsb] depth %d head %.25Rg child %.25Rg cerr %g lerr %g\n",
                             L->depth, head, s, cerr, lerr);
            if (cok) {
                mpfr_add(out, head, s, MPFR_RNDN);
                *err += lerr + cerr;
                ok = true;
            }
        }
        mpfr_clears(head, s, (mpfr_ptr)0);
        mpz_clears(bz, pw, m, (mpz_ptr)0);
        if (ok || x->bad) goto out;
        /* The block series could not be summed (budget): fall through to the
         * honest Levin estimate on the head below. */
    }
    if (L->depth == 0 && shape != NSB_SMOOTH) goto out;
    /* SMOOTH, or a deeper level whose structure could not be followed: Levin on
     * the head.  When the head was cut short the estimate stays, with an error
     * no smaller than the last term times the terms seen. */
    {
        double e = 0.0;
        int want = 0;
        /* Error of the partial sums: the terms' own (block quadrature) error,
         * or the working precision for exact terms. */
        double dS = lerr + nsb_absd(S) * exp2(-(double)x->p + 4.0);
        bool lok = n >= 8 && nsb_levin(x, a, n, k0, dS, s0, out, &e, &want);
        /* Confirm a latest-start estimate whose error is only its predecessor's
         * (twice at most, within the head capacity and the budget). */
        for (int tries = 0; lok && tries < 2 && want > n && want <= cap && !x->budget_hit
                            && e > nsb_absd(out) * exp2(-(double)(x->tb + 2)); tries++) {
            bool st2 = false;
            double t2 = 0.0;
            int n0 = n;
            nsb_fetch(x, L, a, &n, want, k0, S, &lerr, false, &st2, &t2);
            if (x->bad) goto out;
            if (n == n0) break;
            e = 0.0;
            dS = lerr + nsb_absd(S) * exp2(-(double)x->p + 4.0);
            lok = nsb_levin(x, a, n, k0, dS, s0, out, &e, &want);
        }
        if (lok) {
            if (shape != NSB_SMOOTH) {
                double crude = fmax(nsb_absd(a[n - 1]) * (double)n, nsb_powtail(a, n, k0));
                mpfr_t d; mpfr_init2(d, x->p);
                mpfr_sub(d, out, S, MPFR_RNDN);
                double jump = nsb_absd(d);
                mpfr_clear(d);
                if (crude > e) e = crude;
                if (jump > e) e = jump;
            }
            *err += lerr + e;
            ok = true;
            if (nsb_debug())
                mpfr_fprintf(stderr, "[nsb] depth %d levin %.25Rg err %g (lerr %g)\n",
                             L->depth, out, e, lerr);
        } else if (n >= 1 && L->depth > 0) {
            /* Too few terms for Levin (the index budget ends the level).  A
             * geometric tail, when the last two ratios agree on one: the tail
             * a r/(1-r), charged in full plus the ratios' disagreement.  Any
             * other shape: the power-law bound of nsb_powtail. */
            mpfr_set(out, S, MPFR_RNDN);
            double e = nsb_powtail(a, n, k0);
            if (n >= 3 && mpfr_sgn(a[n - 1]) == mpfr_sgn(a[n - 2])
                && mpfr_sgn(a[n - 2]) == mpfr_sgn(a[n - 3])) {
                double r1 = exp2(nsb_log2(a[n - 1]) - nsb_log2(a[n - 2]));
                double r2 = exp2(nsb_log2(a[n - 2]) - nsb_log2(a[n - 3]));
                if (r1 < 0.5 && r2 < 0.5) {
                    mpfr_t t;
                    mpfr_init2(t, x->p);
                    mpfr_mul_d(t, a[n - 1], r1 / (1.0 - r1), MPFR_RNDN);
                    mpfr_add(out, out, t, MPFR_RNDN);
                    e = nsb_absd(t) + nsb_absd(a[n - 1]) * fabs(r1 - r2) / ((1.0 - r1) * (1.0 - r2));
                    mpfr_clear(t);
                }
            } else if (n == 2 && mpfr_sgn(a[1]) == mpfr_sgn(a[0])) {
                /* One ratio only: unconfirmed, so only a steep one is used, and
                 * the tail is charged four times over. */
                double r1 = exp2(nsb_log2(a[1]) - nsb_log2(a[0]));
                if (r1 < 0.1) {
                    mpfr_t t;
                    mpfr_init2(t, x->p);
                    mpfr_mul_d(t, a[1], r1 / (1.0 - r1), MPFR_RNDN);
                    mpfr_add(out, out, t, MPFR_RNDN);
                    e = 4.0 * nsb_absd(t);
                    mpfr_clear(t);
                }
            }
            *err += lerr + e;
            ok = true;
            if (nsb_debug())
                mpfr_fprintf(stderr, "[nsb] depth %d short head n %d: %.25Rg err %g\n",
                             L->depth, n, out, e);
        }
    }
out:
    for (int i = 0; i < n; i++) mpfr_clear(a[i]);
    free(a);
    mpfr_clear(S);
    x->gq = gsave;
    return ok;
}

/* ------------------------------------------------------------------ *
 *  Entry point                                                        *
 * ------------------------------------------------------------------ */

bool nsb_sum(NsbTermFn f, void* ud, const mpz_t n0, const NsbParams* prm,
             mpfr_t out, double* err, int* levels) {
    if (!mpz_fits_slong_p(n0)) return false;
    long k0 = mpz_get_si(n0);
    if (k0 < 1 || k0 > 1000000000L) return false;
    NsbCtx x;
    memset(&x, 0, sizeof x);
    x.f = f; x.ud = ud;
    x.tb = prm->target_bits > 0 ? prm->target_bits : 53;
    x.tq = x.tb + nsb_levin_len(&x) + 16;
    /* The working precision carries tq plus the equispaced interpolant's
     * worst-case Lebesgue constant (~2^q at a one-sided stencil). */
    {
        long need = x.tq + NSB_QMAX + 24;
        x.p = (mpfr_prec_t)(prm->work_bits > need ? prm->work_bits : need);
    }
    x.t_end = mth_mono_seconds() + (prm->budget_s > 0 ? prm->budget_s : 2.0);
    /* Largest index the summand is asked for: 60000 bits (an n of 18000
     * digits), where one evaluation of a digit-length summand still costs
     * microseconds; past it a block level is closed by its tail bound. */
    x.max_index_bits = 60000;
    for (int r = 0; r < NSB_GL_RULES; r++) x.gl_n[r] = 4 << r;
    nsb_bern_build(&x);
    bool ok = nsb_stir_build(&x);
    double e = 0.0;
    if (ok) {
        NsbLevel L0;
        L0.depth = 0;
        L0.parent = NULL;
        mpz_init(L0.b);
        L0.c = 0;
        L0.lmax = (double)x.max_index_bits;
        L0.tq = x.tq;
        x.gq = x.tq;
        mpfr_t s;
        mpfr_init2(s, x.p);
        ok = nsb_sum_level(&x, &L0, k0, s, &e) && !x.bad;
        if (ok) mpfr_set(out, s, MPFR_RNDN);
        mpfr_clear(s);
        mpz_clear(L0.b);
    }
    if (ok) {
        *err = e;
        if (levels) *levels = x.max_level;
    }
    if (nsb_debug())
        fprintf(stderr, "[nsb] %s: %ld summand evaluations, levels %d, p %ld, tq %ld\n",
                ok ? "done" : "declined", x.nevals, x.max_level, (long)x.p, x.tq);
    for (int r = 0; r < NSB_GL_RULES; r++) {
        if (!x.gl_x[r]) continue;
        for (int i = 0; i < x.gl_n[r]; i++) { mpfr_clear(x.gl_x[r][i]); mpfr_clear(x.gl_w[r][i]); }
        free(x.gl_x[r]); free(x.gl_w[r]);
    }
    for (int j = 0; j <= NSB_JMAX; j++) mpfr_clear(x.bern[j]);
    if (x.stir) {
        for (size_t i = 0; i < (size_t)NSB_QE_MAX * NSB_QE_MAX; i++) mpz_clear(x.stir[i]);
        free(x.stir);
    }
    return ok;
}

#endif /* USE_MPFR */

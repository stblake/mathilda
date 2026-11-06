/*
 * nroots_adaptive.c — roots of EXACT integer polynomials with a per-root,
 * adaptively chosen evaluation precision (real Laguerre walk, then Aberth).
 *
 * Why a separate engine.  The fixed-precision engine (nroots_aberth.c) runs
 * every evaluation at one working precision wp = target + guard.  For a
 * polynomial whose monomial expansion cancels heavily on the root set — the
 * classical orthogonal polynomials at high degree, where
 *     sum_k |c_k| |x|^k  /  |p'(x) x|
 * reaches 2^1500 for LaguerreL[1000, x] — p(z) computed at wp is pure rounding
 * noise near most roots, the iteration never meets its tolerance, and the
 * engine spends its whole MaxIterations budget (O(n^2) work per sweep) before
 * returning noise.  Raising wp for every root would be correct but wasteful:
 * the cancellation depends on where the root is (26 digits at x ~ 1, ~470 at
 * x ~ 1000 for L_1000), and most sweeps happen while the iterates are still far
 * from the roots, where low precision is enough.
 *
 * Method.  Every root carries its own evaluation precision; an evaluation
 * also accumulates the rounding-error bound  sum |c_k| |z|^k, and a root is
 * accepted only when its last step is below 2^-(T+6) |z| AND the error bound
 * puts the true root within 2^-(T+2) |z| of it — its T leading bits are
 * justified by the precision actually used.
 *
 *   0. Real walk (ad_real_walk).  Real roots largest first, by Laguerre's
 *      method on p with the found roots divided out implicitly, each step in
 *      MPFR at the precision the current residual needs.  For a real-rooted
 *      polynomial (orthogonal polynomials, Wilkinson-type products) it is
 *      globally and monotonically convergent, so the search for the next root
 *      starts from a safe Newton bound below the last one (or the
 *      extrapolated root sequence) and needs ~3 evaluations per root.  The
 *      walk stops at the first sign of a non-real root.
 *   1. DPE phase, for the roots the walk did not find.  Bini's Newton-polygon
 *      starting points, Aberth with p, p' in double-mantissa / long-exponent
 *      ("DPE") arithmetic, so no coefficient range overflows; a root whose
 *      |p(z)| falls into the rounding-error bound is FROZEN (it cannot be
 *      improved in double).
 *   2. Multiprecision Aberth.  p(z), p'(z) in MPFR at the root's precision by
 *      the real-coefficient Clenshaw (quadratic-divisor) recurrence — 4 real
 *      products per term instead of 8 for complex Horner.  A root inside its
 *      pseudo-zero set has its precision doubled; otherwise the precision is
 *      set for ~24 correct bits of the predicted next residual.  The Aberth
 *      correction only needs a few correct digits (it is a correction), so it
 *      is formed in double from the MPFR values; the repulsion sum uses double
 *      copies of the iterates (MPFR for a pair closer than 2^-40 relatively),
 *      and the walk's roots stay fixed in it.
 *
 * Cost.  A root costs a few MPFR evaluations at roughly the precision its own
 * cancellation demands, O(n M(P_i)) each; roots that need little precision
 * finish early and cheaply.  NSolve[LaguerreL[1000, x] == 0, x,
 * WorkingPrecision -> 50]: ~4000 evaluations at a mean 1400 bits, a few
 * seconds (the fixed-precision engine did not finish in 600 s).
 *
 * Declines (returns -1, the caller falls back to the fixed-precision engine)
 * when a root-radius estimate leaves the double exponent range.
 *
 * Memory: every mpfr_init pairs with a clear before return; the caller owns
 * `roots`.
 */

#include "nroots_internal.h"

#ifdef USE_MPFR

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ *
 *  DPE (double mantissa, long exponent) helpers
 * ------------------------------------------------------------------ */

/* A complex value (re + i im) * 2^e.  Normalised lazily: the mantissas are
 * kept within [2^-400, 2^400] in magnitude (or are zero). */
typedef struct { double re, im; long e; } Dpe;

#define DPE_HI 0x1p400
#define DPE_LO 0x1p-400

static void dpe_norm(Dpe* a) {
    double m = fabs(a->re) > fabs(a->im) ? fabs(a->re) : fabs(a->im);
    if (m == 0.0) { a->e = 0; return; }
    if (m > DPE_HI || m < DPE_LO) {
        int k; (void)frexp(m, &k);
        a->re = ldexp(a->re, -k);
        a->im = ldexp(a->im, -k);
        a->e += k;
    }
}

/* a += b (both DPE). */
static void dpe_add(Dpe* a, const Dpe* b) {
    if (b->re == 0.0 && b->im == 0.0) return;
    if (a->re == 0.0 && a->im == 0.0) { *a = *b; return; }
    long d = b->e - a->e;
    if (d > 0) {
        if (d > 1100) { *a = *b; return; }
        a->re = ldexp(a->re, (int)-d) + b->re;
        a->im = ldexp(a->im, (int)-d) + b->im;
        a->e = b->e;
    } else {
        if (d < -1100) return;
        a->re += ldexp(b->re, (int)d);
        a->im += ldexp(b->im, (int)d);
    }
    dpe_norm(a);
}

/* log2 |a| (or -inf for zero). */
static double dpe_log2abs(const Dpe* a) {
    double m = hypot(a->re, a->im);
    if (m == 0.0) return -INFINITY;
    return log2(m) + (double)a->e;
}

/* ------------------------------------------------------------------ *
 *  Engine state
 * ------------------------------------------------------------------ */

/* Coefficient vector rounded to one evaluation precision (cached by level). */
typedef struct {
    mpfr_prec_t prec;
    mpfr_t*     c;          /* c[0..n] at `prec`, or NULL if not built yet */
} CoefLevel;

typedef struct {
    int          n;
    const mpz_t* z;         /* exact integer coefficients c[0..n]          */
    double*      cm;        /* DPE mantissas of c_k (signed)               */
    long*        ce;        /* DPE exponents                               */
    double*      acm;       /* |cm|                                        */
    CoefLevel*   lev;       /* cache, index = prec / 64                    */
    int          nlev;
    double       log2n2;    /* log2(8 n^2): Clenshaw error-bound factor    */
} AdState;

static mpfr_t* ad_coeffs(AdState* st, mpfr_prec_t P) {
    int idx = (int)(P / 64);
    if (idx >= st->nlev) {
        int nn = idx + 8;
        st->lev = (CoefLevel*)realloc(st->lev, sizeof(CoefLevel) * (size_t)nn);
        for (int i = st->nlev; i < nn; i++) { st->lev[i].c = NULL; st->lev[i].prec = 0; }
        st->nlev = nn;
    }
    CoefLevel* L = &st->lev[idx];
    if (!L->c) {
        L->prec = P;
        L->c = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)(st->n + 1));
        for (int k = 0; k <= st->n; k++) {
            mpfr_init2(L->c[k], P);
            mpfr_set_z(L->c[k], st->z[k], MPFR_RNDN);
        }
    }
    return L->c;
}

/* log2 of  sum_k |c_k| r^k  (r = |z| > 0), by DPE Horner on |c_k|. */
static double ad_log2_absbound(const AdState* st, double r) {
    Dpe a = { 0.0, 0.0, 0 };
    for (int k = st->n; k >= 0; k--) {
        a.re *= r;
        Dpe c = { st->acm[k], 0.0, st->ce[k] };
        dpe_add(&a, &c);
        dpe_norm(&a);
    }
    return dpe_log2abs(&a);
}

/* DPE Horner: p(z), p'(z) for z = zr + i zi (double). */
static void ad_eval_dpe(const AdState* st, double zr, double zi, Dpe* p, Dpe* d) {
    Dpe b = { 0.0, 0.0, 0 }, q = { 0.0, 0.0, 0 };
    for (int k = st->n; k >= 0; k--) {
        /* q = q*z + b */
        double qr = q.re * zr - q.im * zi, qi = q.re * zi + q.im * zr;
        q.re = qr; q.im = qi;
        dpe_add(&q, &b);
        dpe_norm(&q);
        /* b = b*z + c_k */
        double br = b.re * zr - b.im * zi, bi = b.re * zi + b.im * zr;
        b.re = br; b.im = bi;
        Dpe c = { st->cm[k], 0.0, st->ce[k] };
        dpe_add(&b, &c);
        dpe_norm(&b);
    }
    *p = b; *d = q;
}

/* MPFR Clenshaw evaluation of p(z) and p'(z) at z = x + i y, precision P.
 * Real coefficients: with t = 2x, s = x^2 + y^2 (so z^2 = t z - s),
 *     b_k = c_k + t b_{k+1} - s b_{k+2},     p(z)  = b_0 - conj(z) b_1,
 * and the quotient Q of p by (X^2 - t X + s) has coefficients b_2..b_n, so
 *     p'(z) = (z - conj z) Q(z) + b_1,  Q(z) by the same recurrence.
 * Outputs are written to pr, pi, dr, di (initialised by the caller at P). */
static void ad_eval_mp(AdState* st, mpfr_srcptr x, mpfr_srcptr y, mpfr_prec_t P,
                       mpfr_t pr, mpfr_t pi, mpfr_t dr, mpfr_t di) {
    mpfr_t* c = ad_coeffs(st, P);
    int n = st->n;
    mpfr_t t, s, b0, b1, b2, e0, e1, e2, u;
    mpfr_inits2(P, t, s, b0, b1, b2, e0, e1, e2, u, (mpfr_ptr)0);
    mpfr_mul_2ui(t, x, 1, MPFR_RNDN);
    mpfr_sqr(s, x, MPFR_RNDN);
    mpfr_sqr(u, y, MPFR_RNDN);
    mpfr_add(s, s, u, MPFR_RNDN);
    mpfr_set_zero(b1, 1); mpfr_set_zero(b2, 1);
    mpfr_set_zero(e1, 1); mpfr_set_zero(e2, 1);
    for (int k = n; k >= 0; k--) {
        /* b0 = c_k + t*b1 - s*b2 */
        mpfr_mul(b0, t, b1, MPFR_RNDN);
        mpfr_mul(u, s, b2, MPFR_RNDN);
        mpfr_sub(b0, b0, u, MPFR_RNDN);
        mpfr_add(b0, b0, c[k], MPFR_RNDN);
        if (k >= 2) {
            /* Q coefficient q_{k-2} = b_k:  e0 = b_k + t*e1 - s*e2 */
            mpfr_mul(e0, t, e1, MPFR_RNDN);
            mpfr_mul(u, s, e2, MPFR_RNDN);
            mpfr_sub(e0, e0, u, MPFR_RNDN);
            mpfr_add(e0, e0, b0, MPFR_RNDN);
            mpfr_swap(e2, e1); mpfr_swap(e1, e0);
        }
        if (k > 0) { mpfr_swap(b2, b1); mpfr_swap(b1, b0); }
    }
    /* Now b0 = b_0, b1 = b_1, e1 = e_0, e2 = e_1 (Q's Clenshaw terms). */
    /* p(z) = (b0 - x b1) + i (y b1) */
    mpfr_mul(u, x, b1, MPFR_RNDN);
    mpfr_sub(pr, b0, u, MPFR_RNDN);
    mpfr_mul(pi, y, b1, MPFR_RNDN);
    /* Q(z) = (e_0 - x e_1) + i (y e_1);  p' = (b1 - 2y Qi) + i (2y Qr) */
    if (n >= 2) {
        mpfr_t qr, qi;
        mpfr_inits2(P, qr, qi, (mpfr_ptr)0);
        mpfr_mul(u, x, e2, MPFR_RNDN);
        mpfr_sub(qr, e1, u, MPFR_RNDN);
        mpfr_mul(qi, y, e2, MPFR_RNDN);
        mpfr_mul(u, y, qi, MPFR_RNDN);
        mpfr_mul_2ui(u, u, 1, MPFR_RNDN);
        mpfr_sub(dr, b1, u, MPFR_RNDN);
        mpfr_mul(di, y, qr, MPFR_RNDN);
        mpfr_mul_2ui(di, di, 1, MPFR_RNDN);
        mpfr_clears(qr, qi, (mpfr_ptr)0);
    } else {
        mpfr_set(dr, b1, MPFR_RNDN);
        mpfr_set_zero(di, 1);
    }
    mpfr_clears(t, s, b0, b1, b2, e0, e1, e2, u, (mpfr_ptr)0);
}

/* MPFR complex pair -> DPE. */
static void ad_mp_to_dpe(mpfr_srcptr re, mpfr_srcptr im, Dpe* out) {
    long er = 0, ei = 0;
    double mr = mpfr_zero_p(re) ? 0.0 : mpfr_get_d_2exp(&er, re, MPFR_RNDN);
    double mi = mpfr_zero_p(im) ? 0.0 : mpfr_get_d_2exp(&ei, im, MPFR_RNDN);
    if (mr == 0.0 && mi == 0.0) { out->re = out->im = 0.0; out->e = 0; return; }
    long e = (mr == 0.0) ? ei : (mi == 0.0) ? er : (er > ei ? er : ei);
    out->re = (mr == 0.0 || e - er > 1100) ? 0.0 : ldexp(mr, (int)(er - e));
    out->im = (mi == 0.0 || e - ei > 1100) ? 0.0 : ldexp(mi, (int)(ei - e));
    out->e = e;
}

/* N = p / d as DPE (d != 0). */
static void dpe_div(const Dpe* p, const Dpe* d, Dpe* out) {
    double den = d->re * d->re + d->im * d->im;
    out->re = (p->re * d->re + p->im * d->im) / den;
    out->im = (p->im * d->re - p->re * d->im) / den;
    out->e = p->e - d->e;
    dpe_norm(out);
}

/* ------------------------------------------------------------------ *
 *  Newton-polygon starting points (Bini 1996), on log2|c_k|.
 * ------------------------------------------------------------------ */
static int ad_initial_points(const AdState* st, double* zr, double* zi) {
    int n = st->n;
    double* lk = (double*)malloc(sizeof(double) * (size_t)(n + 1));
    int* kk = (int*)malloc(sizeof(int) * (size_t)(n + 1));
    int* hull = (int*)malloc(sizeof(int) * (size_t)(n + 1));
    int np = 0, hn = 0, ok = 1;
    for (int k = 0; k <= n; k++) {
        if (st->cm[k] == 0.0) continue;
        lk[np] = log2(fabs(st->cm[k])) + (double)st->ce[k];
        kk[np] = k; np++;
    }
    for (int i = 0; i < np; i++) {
        while (hn >= 2) {
            int a = hull[hn - 2], b = hull[hn - 1];
            double cross = (double)(kk[b] - kk[a]) * (lk[i] - lk[a])
                         - (lk[b] - lk[a]) * (double)(kk[i] - kk[a]);
            if (cross >= 0.0) hn--; else break;
        }
        hull[hn++] = i;
    }
    int idx = 0;
    const double sigma = 0.7;
    for (int e = 0; e + 1 < hn; e++) {
        int ia = hull[e], ib = hull[e + 1];
        int dk = kk[ib] - kk[ia];
        if (dk <= 0) continue;
        double l2r = (lk[ia] - lk[ib]) / (double)dk;
        if (l2r > 900.0 || l2r < -900.0) { ok = 0; break; }
        double r = exp2(l2r);
        double ph = 2.0 * M_PI * (double)e / (double)n;
        for (int t = 0; t < dk && idx < n; t++) {
            double th = (2.0 * M_PI / (double)dk) * (double)t + sigma + ph;
            zr[idx] = r * cos(th); zi[idx] = r * sin(th); idx++;
        }
    }
    for (; idx < n; idx++) {
        double th = 2.0 * M_PI * (double)idx / (double)n + sigma;
        zr[idx] = cos(th); zi[idx] = sin(th);
    }
    free(lk); free(kk); free(hull);
    return ok;
}

/* S_i = sum_{j != i} 1 / (z_i - z_j) in double; pairs closer than 2^-40
 * relative use the MPFR iterates (`xr`, `xi`) to avoid cancellation. */
static void ad_repulsion(int n, int i, const double* zr, const double* zi,
                         mpfr_t* xr, mpfr_t* xi, double* sr, double* si) {
    double ar = 0.0, ai = 0.0;
    double zabs = hypot(zr[i], zi[i]);
    for (int j = 0; j < n; j++) {
        if (j == i) continue;
        double dr = zr[i] - zr[j], di = zi[i] - zi[j];
        double m2 = dr * dr + di * di;
        if (xr && m2 < 0x1p-80 * zabs * zabs) {
            /* Close pair: difference in MPFR, then to double. */
            mpfr_t a; mpfr_init2(a, mpfr_get_prec(xr[i]));
            mpfr_sub(a, xr[i], xr[j], MPFR_RNDN); dr = mpfr_get_d(a, MPFR_RNDN);
            mpfr_sub(a, xi[i], xi[j], MPFR_RNDN); di = mpfr_get_d(a, MPFR_RNDN);
            mpfr_clear(a);
            m2 = dr * dr + di * di;
        }
        if (m2 == 0.0 || !isfinite(m2)) continue;
        ar += dr / m2; ai -= di / m2;
    }
    *sr = ar; *si = ai;
}

/* w = N / (1 - N S) in DPE (N DPE, S double). */
static void ad_aberth_corr(const Dpe* N, double sr, double si, Dpe* w) {
    Dpe ns = { -(N->re * sr - N->im * si), -(N->re * si + N->im * sr), N->e };
    dpe_norm(&ns);
    Dpe den = { 1.0, 0.0, 0 };
    dpe_add(&den, &ns);                      /* den = 1 - N S */
    if ((den.re == 0.0 && den.im == 0.0) || !isfinite(den.re) || !isfinite(den.im)) {
        *w = *N; return;
    }
    dpe_div(N, &den, w);
}

/* ------------------------------------------------------------------ *
 *  Real roots in decreasing order: Laguerre with implicit deflation
 * ------------------------------------------------------------------ */

/* Real MPFR Horner of p, p' and (when h != NULL) p''/2 at x, precision P.
 * x carries only the root precision Tp << P, so each product is a short
 * P x Tp multiplication; mul + add (not mpfr_fma, which allocates an exact
 * product per call) keeps the inner loop allocation-free. */
static void ad_eval_real_mp(AdState* st, mpfr_srcptr x, mpfr_prec_t P,
                            mpfr_t p, mpfr_t d, mpfr_t h) {
    mpfr_t* c = ad_coeffs(st, P);
    mpfr_set(p, c[st->n], MPFR_RNDN);
    mpfr_set_zero(d, 1);
    if (h) mpfr_set_zero(h, 1);
    for (int k = st->n - 1; k >= 0; k--) {
        if (h) { mpfr_mul(h, h, x, MPFR_RNDN); mpfr_add(h, h, d, MPFR_RNDN); }
        mpfr_mul(d, d, x, MPFR_RNDN); mpfr_add(d, d, p, MPFR_RNDN);
        mpfr_mul(p, p, x, MPFR_RNDN); mpfr_add(p, p, c[k], MPFR_RNDN);
    }
}

static double ad_mpfr_log2abs(mpfr_srcptr a) {
    if (mpfr_zero_p(a)) return -INFINITY;
    long e = 0;
    double m = mpfr_get_d_2exp(&e, a, MPFR_RNDN);
    return log2(fabs(m)) + (double)e;
}

/* a / b as a double scaled by 2^sc (b != 0); overflow-safe via exponents. */
static double ad_ratio_scaled(mpfr_srcptr a, mpfr_srcptr b, long sc) {
    if (mpfr_zero_p(a)) return 0.0;
    long ea = 0, eb = 0;
    double ma = mpfr_get_d_2exp(&ea, a, MPFR_RNDN);
    double mb = mpfr_get_d_2exp(&eb, b, MPFR_RNDN);
    long e = ea - eb + sc;
    if (e > 1000) e = 1000;
    if (e < -1100) return 0.0;
    return ldexp(ma / mb, (int)e);
}

/* For a real polynomial whose roots are all real (the orthogonal
 * polynomials, characteristic polynomials of symmetric matrices, ...),
 * Aberth's simultaneous iteration from Newton-polygon starting points is the
 * wrong tool when the expansion cancels: the tropical radii are far from the
 * roots, the iterates crawl in along arcs for ~100 sweeps, and every sweep
 * already needs the high precision the cancellation demands.  Laguerre's
 * method on q = p / prod (x - r_j) (the roots found so far divided out
 * implicitly: q'/q = p'/p - sum 1/(x - r_j)) instead converges globally and
 * monotonically for a real-rooted polynomial — from above the largest
 * remaining root it moves down onto it, cubically, in a handful of steps.
 * Started just below each root found, it walks the real roots in decreasing
 * order, each evaluated at the precision its own cancellation requires.
 *
 * Every accepted root is certified (step and evaluation-error uncertainty
 * below the target) and strictly below the previous one, so the found roots
 * are distinct genuine roots.  The walk stops at the first sign of a non-real
 * root (Laguerre's radicand mH - G^2 is >= 0 at every real point exactly when
 * all remaining roots are real), at non-convergence, or below the root bound;
 * the caller completes the remaining roots with Aberth.  Returns the number
 * of roots written to rr[0..] (each mpfr_init2'd by the caller at Tp). */
static int ad_real_walk(AdState* st, long T, mpfr_prec_t Tp, mpfr_t* rr) {
    int n = st->n;
    double tolw = -(double)(T + 6), tolu = -(double)(T + 2);
    double lerr = log2(4.0 * (double)(n + 1));     /* Horner error factor */

    /* Fujiwara bound: every root modulus is below 2 max |c_k/c_n|^(1/(n-k)). */
    double l2cn = log2(fabs(st->cm[n])) + (double)st->ce[n];
    double lB = -INFINITY;
    for (int k = 0; k < n; k++) {
        if (st->cm[k] == 0.0) continue;
        double v = (log2(fabs(st->cm[k])) + (double)st->ce[k] - l2cn) / (double)(n - k);
        if (v > lB) lB = v;
    }
    lB += 1.0;
    if (!(lB < 900.0 && lB > -900.0)) return 0;

    mpfr_t x, p, d, h, a, G, g, hh;
    mpfr_prec_t Pc = Tp + 64;                       /* step arithmetic   */
    mpfr_init2(x, Tp); mpfr_inits2(64, p, d, h, a, G, g, hh, (mpfr_ptr)0);
    mpfr_set_d(x, exp2(lB) * 1.0625, MPFR_RNDU);
    double* rd = (double*)malloc(sizeof(double) * (size_t)n);
    int found = 0;
    mpfr_prec_t P = 64, P_first = 64;

    while (found < n) {
        int m = n - found;
        int ok = 0;
        P_first = 0;
        for (int it = 0; it < 100; it++) {
            mpfr_set_prec(p, P); mpfr_set_prec(d, P); mpfr_set_prec(h, P);
            ad_eval_real_mp(st, x, P, p, d, h);
            double xd = mpfr_get_d(x, MPFR_RNDN);
            double lx = log2(fabs(xd));
            double lab = ad_log2_absbound(st, fabs(xd));
            double leb = lab - (double)P + lerr;
            double lp = ad_mpfr_log2abs(p), ld = ad_mpfr_log2abs(d);
            if (ld <= leb - lx + 8.0) { P *= 2; if (P > 4000000) break; continue; }
            if (lp <= leb + 8.0) {
                if (leb - ld - lx <= tolu) { ok = 1; break; }
                P *= 2; if (P > 4000000) break;
                continue;
            }
            /* Laguerre on q = p / prod (x - r_j), in MPFR at the root
             * precision (a step formed in double would gain only ~50 bits a
             * iteration near the root):
             *   g = p'/p - S1,   hh = (p'/p)^2 - p''/p - S2,
             *   S1 = sum 1/(x - r_j),  S2 = sum 1/(x - r_j)^2  (double: the
             *   found roots are a root gap away, their terms are not dominant),
             *   step = m / (g + sgn(g) sqrt((m-1)(m hh - g^2))). */
            double S1 = 0.0, S2 = 0.0;
            for (int j = 0; j < found; j++) {
                double df = xd - rd[j];
                if (fabs(df) < 0x1p-30 * fabs(xd)) {
                    mpfr_set_prec(a, Tp);
                    mpfr_sub(a, x, rr[j], MPFR_RNDN);
                    df = mpfr_get_d(a, MPFR_RNDN);
                }
                double t = 1.0 / df;
                S1 += t; S2 += t * t;
            }
            mpfr_set_prec(G, Pc); mpfr_set_prec(g, Pc); mpfr_set_prec(hh, Pc);
            mpfr_set_prec(a, Pc);
            mpfr_div(G, d, p, MPFR_RNDN);               /* p'/p           */
            mpfr_sub_d(g, G, S1, MPFR_RNDN);
            /* Dividing the found roots out cancels log2|G/g| bits of p'/p;
             * p'/p must carry those plus ~20 more, else re-evaluate higher. */
            {
                double acc = lp - leb;
                if (ld + lx - leb < acc) acc = ld + lx - leb;
                double canc = mpfr_zero_p(g) ? 64.0
                            : ad_mpfr_log2abs(G) - ad_mpfr_log2abs(g);
                if (canc < 0.0) canc = 0.0;
                if (acc < canc + 20.0) {
                    mpfr_prec_t np = P + (mpfr_prec_t)ceil(canc + 20.0 - acc) + 32;
                    P = (np + 63) / 64 * 64;
                    if (P > 4000000) break;
                    continue;
                }
            }
            mpfr_sqr(hh, G, MPFR_RNDN);
            mpfr_div(a, h, p, MPFR_RNDN);
            mpfr_mul_2ui(a, a, 1, MPFR_RNDN);           /* p''/p          */
            mpfr_sub(hh, hh, a, MPFR_RNDN);
            mpfr_sub_d(hh, hh, S2, MPFR_RNDN);
            /* rad = (m-1)(m hh - g^2) */
            mpfr_mul_ui(hh, hh, (unsigned long)m, MPFR_RNDN);
            mpfr_sqr(a, g, MPFR_RNDN);
            int gsq_exp_ok = !mpfr_zero_p(a);
            mpfr_sub(hh, hh, a, MPFR_RNDN);
            if (mpfr_sgn(hh) < 0) {
                /* Rounding may leave a tiny negative radicand at a simple
                 * root; a substantial one means non-real roots remain. */
                if (gsq_exp_ok) mpfr_mul_2si(a, a, -26, MPFR_RNDN);
                if (mpfr_cmpabs(hh, a) > 0) break;
                mpfr_set_zero(hh, 1);
            }
            mpfr_mul_ui(hh, hh, (unsigned long)(m - 1), MPFR_RNDN);
            mpfr_sqrt(hh, hh, MPFR_RNDN);
            if (mpfr_sgn(g) >= 0) mpfr_add(a, g, hh, MPFR_RNDN);
            else                  mpfr_sub(a, g, hh, MPFR_RNDN);
            if (mpfr_zero_p(a) || !mpfr_number_p(a)) break;
            mpfr_ui_div(a, (unsigned long)m, a, MPFR_RNDN); /* step       */
            if (!mpfr_number_p(a)) break;
            double rel = ad_mpfr_log2abs(a) - lx;
            if (P_first == 0) P_first = P;
            mpfr_sub(x, x, a, MPFR_RNDN);
            if (found > 0 && mpfr_cmp(x, rr[found - 1]) >= 0) break;
            if (log2(fabs(mpfr_get_d(x, MPFR_RNDN))) > lB + 1.0) break;
            double lunc = leb - ld - lx;
            if (rel <= tolw && lunc <= tolu) { ok = 1; break; }
            double rel_next = (rel < -4.0) ? 3.0 * rel : rel;
            if (rel_next < -(double)(T + 4)) rel_next = -(double)(T + 4);
            double lp_next = (rel < -4.0) ? ld + lx + rel_next : lp;
            double needd = lab - lp_next + lerr + ((rel_next <= -(double)(T + 4)) ? 4.0 : 24.0);
            mpfr_prec_t need = (needd < 64.0) ? 64 : (mpfr_prec_t)ceil(needd);
            need = (need + 63) / 64 * 64;
            if (need > P || rel < -4.0) P = need;
        }
        if (!ok) break;
        if (found > 0 && mpfr_cmp(x, rr[found - 1]) >= 0) break;
        mpfr_set(rr[found], x, MPFR_RNDN);
        rd[found] = mpfr_get_d(x, MPFR_RNDN);
        found++;
        if (found == n) break;
        /* Restart point.  The safe one is a Newton step from the root on the
         * newly deflated q = p / prod (x - r_j), j <= found:  at the root
         *     q'/q = p''/(2 p') - sum_{j < found} 1/(r - r_j),
         * and for a real-rooted q (every remaining root below r) Newton from
         * above never overshoots the largest remaining root, so 1/(q'/q) is a
         * LOWER bound on the gap to the next root.  Root sequences of the
         * polynomials this path serves are smooth, so the quadratic
         * extrapolation 3 r_k - 3 r_{k-1} + r_{k-2} is usually far closer;
         * it is used when it lies below the Newton point (consistent with the
         * bound), else the Newton point.  A bad guess only costs iterations:
         * Laguerre converges to the root nearest the start, and a skipped
         * root shows up as a short count, never as a wrong one. */
        {
            double xr0 = rd[found - 1];
            double gq = ad_ratio_scaled(h, d, 0);
            for (int j = 0; j < found - 1; j++) gq -= 1.0 / (xr0 - rd[j]);
            if (!(gq > 0.0) || !isfinite(gq)) break;   /* not real-rooted */
            double xn = xr0 - 1.0 / gq;
            double xs = xn;
            if (found >= 3) {
                double xe = 3.0 * rd[found - 1] - 3.0 * rd[found - 2] + rd[found - 3];
                double gap = rd[found - 2] - rd[found - 1];
                /* Pad the guess upward by the size of the last third
                 * difference: Laguerre from just above the root is monotone. */
                double third = fabs(rd[found - 1] - 3.0 * rd[found - 2]
                                    + 3.0 * rd[found - 3] - (found >= 4 ? rd[found - 4] : rd[found - 3]));
                if (found < 4) third = 0.01 * gap;
                xe += 2.0 * third;
                if (xe < xn && xe > xr0 - 4.0 * gap) xs = xe;
            }
            mpfr_set_prec(a, 64);
            mpfr_set_d(a, xr0 - xs, MPFR_RNDN);
            mpfr_sub(x, x, a, MPFR_RNDN);
        }
        /* Next root's precision: the first precision that sufficed for this
         * one (the cancellation varies slowly along the root set). */
        if (P_first > 0) P = P_first;
    }
    mpfr_clears(x, p, d, h, a, G, g, hh, (mpfr_ptr)0);
    free(rd);
    return found;
}

/* ------------------------------------------------------------------ *
 *  Driver
 * ------------------------------------------------------------------ */
int nr_aberth_exact(const mpz_t* coef, int n, long target_bits,
                    ncpx* roots, mpfr_prec_t wp, int* uncertified) {
    if (uncertified) *uncertified = 0;
    if (n < 1 || mpz_sgn(coef[n]) == 0 || mpz_sgn(coef[0]) == 0) return -1;

    AdState st;
    memset(&st, 0, sizeof st);
    st.n = n; st.z = coef;
    st.cm = (double*)malloc(sizeof(double) * (size_t)(n + 1));
    st.acm = (double*)malloc(sizeof(double) * (size_t)(n + 1));
    st.ce = (long*)malloc(sizeof(long) * (size_t)(n + 1));
    for (int k = 0; k <= n; k++) {
        long e = 0;
        double m = mpz_sgn(coef[k]) ? mpz_get_d_2exp(&e, coef[k]) : 0.0;
        st.cm[k] = m; st.acm[k] = fabs(m); st.ce[k] = e;
    }
    st.log2n2 = log2(8.0 * (double)n * (double)n);

    long T = target_bits;
    double tolu = -(double)(T + 2);      /* log2 relative uncertainty bound  */
    mpfr_prec_t Tp = wp > (mpfr_prec_t)(T + 16) ? wp : (mpfr_prec_t)(T + 16);
    mpfr_t* xr = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)n);
    mpfr_t* xi = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)n);
    for (int i = 0; i < n; i++) {
        mpfr_init2(xr[i], Tp); mpfr_init2(xi[i], Tp);
        mpfr_set_zero(xi[i], 1);
    }

    /* ---- 0. Real roots, largest first, while the polynomial looks
     * real-rooted (all of them for an orthogonal polynomial). ---- */
    int nfound = ad_real_walk(&st, T, Tp, xr);

    double* zr = (double*)malloc(sizeof(double) * (size_t)n);
    double* zi = (double*)malloc(sizeof(double) * (size_t)n);
    int* done = (int*)calloc((size_t)n, sizeof(int));
    int* frozen = (int*)calloc((size_t)n, sizeof(int));
    for (int i = 0; i < nfound; i++) {
        done[i] = 1;
        zr[i] = mpfr_get_d(xr[i], MPFR_RNDN); zi[i] = 0.0;
    }
    if (nfound < n) {
        /* Aberth for the rest, the found roots held fixed (they repel the
         * free iterates through S_i).  Starting points: the first n - nfound
         * of the Newton-polygon placement. */
        double* ir = (double*)malloc(sizeof(double) * (size_t)n);
        double* ii = (double*)malloc(sizeof(double) * (size_t)n);
        int okp = ad_initial_points(&st, ir, ii);
        for (int i = nfound; i < n; i++) { zr[i] = ir[i - nfound]; zi[i] = ii[i - nfound]; }
        free(ir); free(ii);
        if (!okp) {
            for (int i = 0; i < n; i++) { mpfr_clear(xr[i]); mpfr_clear(xi[i]); }
            free(xr); free(xi); free(zr); free(zi); free(done); free(frozen);
            for (int l = 0; l < st.nlev; l++) {
                if (!st.lev[l].c) continue;
                for (int k = 0; k <= n; k++) mpfr_clear(st.lev[l].c[k]);
                free(st.lev[l].c);
            }
            free(st.lev);
            free(st.cm); free(st.acm); free(st.ce);
            return -1;
        }
    }
    double tolw = -(double)(T + 6);      /* log2 relative correction bound   */
    /* DPE freeze threshold: p(z) must carry at least ~10 significant bits
     * above the (typical, sqrt(n)-scaled) rounding error of the double Horner
     * sum for the correction to mean anything. */
    double u_dpe = 0.5 * log2((double)(n + 1)) - 53.0 + 10.0;

    /* ---- 1. DPE phase. ---- */
    int max_dpe = 60 + n / 4;
    for (int sweep = 0; sweep < max_dpe; sweep++) {
        int active = 0;
        for (int i = 0; i < n; i++) {
            if (done[i] || frozen[i]) continue;
            Dpe p, d, N, w;
            ad_eval_dpe(&st, zr[i], zi[i], &p, &d);
            double za = hypot(zr[i], zi[i]);
            double lab = ad_log2_absbound(&st, za);
            double lp = dpe_log2abs(&p), ld = dpe_log2abs(&d);
            double leb = lab + u_dpe;
            if (lp <= leb + 1.0 || ld == -INFINITY) { frozen[i] = 1; continue; }
            active++;
            dpe_div(&p, &d, &N);
            double sr, si;
            ad_repulsion(n, i, zr, zi, NULL, NULL, &sr, &si);
            ad_aberth_corr(&N, sr, si, &w);
            double wr = ldexp(w.re, (int)(w.e > 1100 ? 1100 : w.e < -1100 ? -1100 : w.e));
            double wi = ldexp(w.im, (int)(w.e > 1100 ? 1100 : w.e < -1100 ? -1100 : w.e));
            if (!isfinite(wr) || !isfinite(wi)) { frozen[i] = 1; continue; }
            zr[i] -= wr; zi[i] -= wi;
            double lw = dpe_log2abs(&w);
            double lz = log2(hypot(zr[i], zi[i]));
            if (lw - lz <= -40.0) frozen[i] = 1;   /* hand over to MPFR */
        }
        if (active == 0) break;
    }

    /* ---- 2. Multiprecision phase. ---- */
    mpfr_prec_t* P = (mpfr_prec_t*)malloc(sizeof(mpfr_prec_t) * (size_t)n);
    int* iters = (int*)calloc((size_t)n, sizeof(int));
    for (int i = 0; i < n; i++) {
        if (!done[i]) {
            mpfr_set_d(xr[i], zr[i], MPFR_RNDN);
            mpfr_set_d(xi[i], zi[i], MPFR_RNDN);
        }
        P[i] = 64;
    }
    /* Cap on the evaluation precision: far beyond any cancellation an integer
     * polynomial of this size can exhibit (Mahler-type separation bounds are
     * O(n (log n + coefficient bits))). */
    double maxbits = 0.0;
    for (int k = 0; k <= n; k++) {
        double b = (st.cm[k] != 0.0) ? (double)st.ce[k] : 0.0;
        if (b > maxbits) maxbits = b;
    }
    double pcap_d = 4.0 * ((double)n * (log2((double)n) + 2.0) + 2.0 * maxbits) + (double)T + 256.0;
    if (pcap_d > 4.0e6) pcap_d = 4.0e6;
    mpfr_prec_t Pcap = (mpfr_prec_t)pcap_d;
    int max_it = 80;

    mpfr_t pr, pi, dr, di;
    mpfr_inits2(64, pr, pi, dr, di, (mpfr_ptr)0);
    int remaining = 0;
    for (int i = 0; i < n; i++) if (!done[i]) remaining++;
    while (remaining > 0) {
        int progressed = 0;
        for (int i = 0; i < n; i++) {
            if (done[i]) continue;
            if (iters[i]++ > max_it || P[i] > Pcap) {
                done[i] = 1; remaining--;
                if (uncertified) (*uncertified)++;
                continue;
            }
            progressed = 1;
            mpfr_prec_t Pi = P[i];
            mpfr_set_prec(pr, Pi); mpfr_set_prec(pi, Pi);
            mpfr_set_prec(dr, Pi); mpfr_set_prec(di, Pi);
            ad_eval_mp(&st, xr[i], xi[i], Pi, pr, pi, dr, di);
            Dpe p, d, N, w;
            ad_mp_to_dpe(pr, pi, &p);
            ad_mp_to_dpe(dr, di, &d);
            double za = hypot(zr[i], zi[i]);
            double lz = log2(za);
            double lab = ad_log2_absbound(&st, za);
            double leb = lab - (double)Pi + st.log2n2;
            double lp = dpe_log2abs(&p), ld = dpe_log2abs(&d);
            if (ld == -INFINITY || ld <= leb - lz + 8.0) {
                /* p' unresolved: double the precision. */
                P[i] = Pi * 2; continue;
            }
            if (lp <= leb + 8.0) {
                /* Inside (or at the edge of) the pseudo-zero set at P_i: the
                 * root is located to within eb/|p'|.  Accept if that meets
                 * the target, else raise the precision. */
                if (leb - ld - lz <= tolu) { done[i] = 1; remaining--; continue; }
                P[i] = Pi * 2; continue;
            }
            dpe_div(&p, &d, &N);
            double sr, si;
            ad_repulsion(n, i, zr, zi, xr, xi, &sr, &si);
            ad_aberth_corr(&N, sr, si, &w);
            /* x -= w */
            {
                mpfr_t a; mpfr_init2(a, 64);
                if (w.re != 0.0) {
                    mpfr_set_d(a, w.re, MPFR_RNDN); mpfr_mul_2si(a, a, w.e, MPFR_RNDN);
                    mpfr_sub(xr[i], xr[i], a, MPFR_RNDN);
                }
                if (w.im != 0.0) {
                    mpfr_set_d(a, w.im, MPFR_RNDN); mpfr_mul_2si(a, a, w.e, MPFR_RNDN);
                    mpfr_sub(xi[i], xi[i], a, MPFR_RNDN);
                }
                mpfr_clear(a);
            }
            zr[i] = mpfr_get_d(xr[i], MPFR_RNDN);
            zi[i] = mpfr_get_d(xi[i], MPFR_RNDN);
            double lw = dpe_log2abs(&w);
            double rel = lw - lz;
            double lunc = leb - ld - lz;         /* log2 relative uncertainty */
            if (rel <= tolw && lunc <= tolu) { done[i] = 1; remaining--; continue; }
            /* Precision for the next step: enough for ~24 correct bits of the
             * predicted next residual |p| ~ |p'| |z| 2^rel' (cubic convergence
             * once the correction is small), never more than the target needs. */
            double rel_next = (rel < -4.0) ? 3.0 * rel : rel;
            if (rel_next < -(double)(T + 8)) rel_next = -(double)(T + 8);
            double lp_next = (rel < -4.0) ? ld + lz + rel_next : lp;
            double needd = lab - lp_next + st.log2n2 + 24.0;
            mpfr_prec_t need = (needd < 64.0) ? 64 : (mpfr_prec_t)ceil(needd);
            need = (need + 63) / 64 * 64;
            if (need > Pi || rel < -4.0) P[i] = need;
        }
        if (!progressed) break;
    }
    mpfr_clears(pr, pi, dr, di, (mpfr_ptr)0);

    for (int i = 0; i < n; i++) {
        mpfr_set(roots[i].re, xr[i], MPFR_RNDN);
        mpfr_set(roots[i].im, xi[i], MPFR_RNDN);
        mpfr_clear(xr[i]); mpfr_clear(xi[i]);
    }
    free(xr); free(xi); free(P); free(iters);
    for (int l = 0; l < st.nlev; l++) {
        if (!st.lev[l].c) continue;
        for (int k = 0; k <= n; k++) mpfr_clear(st.lev[l].c[k]);
        free(st.lev[l].c);
    }
    free(st.lev);
    free(zr); free(zi); free(done); free(frozen);
    free(st.cm); free(st.acm); free(st.ce);
    return 0;
}

#endif /* USE_MPFR */

/* sf_machine_meijerg.c -- machine-precision Meijer G-function.
 *
 * See sf_machine_meijerg.h for the method (Slater's theorem with regularized
 * pFq terms) and the decline contract. The per-term pFq values come from the
 * shared hypergeometric kernels (sf_machine_pfq / sf_machine_pfq_c), so every
 * Slater term is evaluated by exactly the code HypergeometricPFQ uses.
 */
#include "sf_machine_meijerg.h"
#include "sf_machine.h"        /* sf_machine_pfq (real)              */
#include "sf_machine_hyper.h"  /* sf_machine_pfq_c (complex)         */
#include "sf_machine_gamma.h"  /* sfg_gamma_c, sfg_loggamma_c        */
#include "flint_num_bridge.h"  /* flint_mg_eval_d (arb fallback), flint_mg_laurent */
#include "sf_dd.h"             /* double-double: the log-case table and re-sum */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef MATHILDA_THREADS
#include <pthread.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2 0.693147180559945309417
#endif

/* ------------------------------------------------------------------ */
/*  (re, im) arithmetic                                                */
/* ------------------------------------------------------------------ */

static inline void c_mul(double ar, double ai, double br, double bi,
                         double* rr, double* ri) {
    double r = ar * br - ai * bi;
    double i = ar * bi + ai * br;
    *rr = r; *ri = i;
}

static inline void c_div(double ar, double ai, double br, double bi,
                         double* rr, double* ri) {
    /* Smith's algorithm: no overflow for large |b|. */
    if (fabs(br) >= fabs(bi)) {
        double t = bi / br, d = br + bi * t;
        *rr = (ar + ai * t) / d;
        *ri = (ai - ar * t) / d;
    } else {
        double t = br / bi, d = br * t + bi;
        *rr = (ar * t + ai) / d;
        *ri = (ai * t - ar) / d;
    }
}

static inline double c_abs(double r, double i) { return hypot(r, i); }

static inline void c_exp(double zr, double zi, double* rr, double* ri) {
    double m = exp(zr);
    if (zi == 0.0) { *rr = m; *ri = 0.0; return; }
    *rr = m * cos(zi);
    *ri = m * sin(zi);
}

/* Principal log; a signed-zero imaginary part is read as +0, so the negative
 * real axis has arg +pi (Mathematica's branch, continuous from above). */
static inline void c_log(double zr, double zi, double* lr, double* li) {
    if (zi == 0.0) zi = 0.0;               /* -0.0 -> +0.0 */
    *lr = log(c_abs(zr, zi));
    *li = atan2(zi, zr);
}

/* sin(pi x) and cos(pi x) with the argument reduced EXACTLY to |f| <= 1/4 first,
 * so sin(pi k) is exactly 0 and the pi/sin(pi d) factors are as accurate as the
 * difference d itself. */
static void sincospi(double x, double* s, double* c) {
    double n = nearbyint(2.0 * x);
    double f = x - 0.5 * n;                /* exact: |f| <= 1/4 */
    double sf = sin(M_PI * f), cf = cos(M_PI * f);
    double k4 = fmod(n, 4.0);
    if (k4 < 0) k4 += 4.0;
    switch ((int)k4) {
        case 0:  *s = sf;  *c = cf;  break;
        case 1:  *s = cf;  *c = -sf; break;
        case 2:  *s = -sf; *c = -cf; break;
        default: *s = -cf; *c = sf;  break;
    }
}

/* sin(pi (x + i y)). */
static void c_sinpi(double x, double y, double* rr, double* ri) {
    double s, c;
    sincospi(x, &s, &c);
    if (y == 0.0) { *rr = s; *ri = 0.0; return; }
    *rr = s * cosh(M_PI * y);
    *ri = c * sinh(M_PI * y);
}

static inline bool is_nonpos_int(double r, double i) {
    return i == 0.0 && r <= 0.0 && r == nearbyint(r);
}

static inline bool is_int(double r, double i) {
    return i == 0.0 && r == nearbyint(r);
}

/* ------------------------------------------------------------------ */
/*  The coefficient of one Slater term                                 */
/* ------------------------------------------------------------------ */

/* A running product kept as a finite complex mantissa plus a power-of-two
 * exponent, so a coefficient that is the ratio of two huge Gamma values stays
 * representable until the very end. */
typedef struct { double r, i; long e; bool zero; } Prod;

static void prod_init(Prod* P) { P->r = 1.0; P->i = 0.0; P->e = 0; P->zero = false; }

static void prod_norm(Prod* P) {
    double a = c_abs(P->r, P->i);
    if (a == 0.0 || !isfinite(a)) return;
    int ex;
    (void)frexp(a, &ex);
    if (ex > 64 || ex < -64) {
        P->r = ldexp(P->r, -ex); P->i = ldexp(P->i, -ex); P->e += ex;
    }
}

static void prod_mul(Prod* P, double r, double i) {
    c_mul(P->r, P->i, r, i, &P->r, &P->i);
    prod_norm(P);
}

static void prod_div(Prod* P, double r, double i) {
    c_div(P->r, P->i, r, i, &P->r, &P->i);
    prod_norm(P);
}

/* Multiply by Gamma(z) (inv = false) or 1/Gamma(z) (inv = true). A pole of the
 * numerator Gamma is an inconsistent parameter set: returns false. A pole under
 * 1/Gamma zeroes the product. Large arguments, where Gamma itself overflows,
 * go through LogGamma. */
static bool prod_gamma(Prod* P, double zr, double zi, bool inv) {
    if (is_nonpos_int(zr, zi)) {
        if (inv) { P->zero = true; return true; }
        return false;
    }
    double gr, gi;
    if (c_abs(zr, zi) < 150.0 && sfg_gamma_c(zr, zi, &gr, &gi) &&
        isfinite(gr) && isfinite(gi) && (gr != 0.0 || gi != 0.0)) {
        if (inv) prod_div(P, gr, gi); else prod_mul(P, gr, gi);
        return true;
    }
    double lr, li;
    if (!sfg_loggamma_c(zr, zi, &lr, &li) || !isfinite(lr) || !isfinite(li)) return false;
    if (inv) { lr = -lr; li = -li; }
    /* exp(lr + i li) = 2^(lr / ln 2) * e^(i li): split off the binary exponent. */
    double l2 = lr / M_LN2;
    double ie = floor(l2);
    double mag = exp((l2 - ie) * M_LN2);
    double er = mag * cos(li), ei = mag * sin(li);
    c_mul(P->r, P->i, er, ei, &P->r, &P->i);
    P->e += (long)ie;
    prod_norm(P);
    return true;
}

static bool prod_value(const Prod* P, double* rr, double* ri) {
    if (P->zero) { *rr = 0.0; *ri = 0.0; return true; }
    if (P->e > 2000 || P->e < -2000) return false;
    *rr = ldexp(P->r, (int)P->e);
    *ri = ldexp(P->i, (int)P->e);
    return isfinite(*rr) && isfinite(*ri);
}

/* ------------------------------------------------------------------ */
/*  One Slater expansion                                               */
/* ------------------------------------------------------------------ */

/* Build S for G^{m,n}_{p,q}(A; B; .) expanded about 0 (the "z-series").
 * Returns false when the parameters are inconsistent (a numerator Gamma pole:
 * a_j - b_h a positive integer for j <= n, h <= m). */
static bool build_series(SfMgSeries* S,
                         const double* Ar, const double* Ai, size_t n, size_t p,
                         const double* Br, const double* Bi, size_t m, size_t q) {
    memset(S, 0, sizeof *S);
    S->nt = m;
    S->pu = p;
    S->ql = q ? q - 1 : 0;
    S->sgn = (((long)p - (long)m - (long)n) % 2 != 0) ? -1.0 : 1.0;
    S->real = true;
    for (size_t j = 0; j < p; j++) if (Ai[j] != 0.0) S->real = false;
    for (size_t j = 0; j < q; j++) if (Bi[j] != 0.0) S->real = false;

    /* The log case: two of b_1..b_m an integer apart. */
    for (size_t h = 0; h < m; h++)
        for (size_t j = h + 1; j < m; j++)
            if (is_int(Br[j] - Br[h], Bi[j] - Bi[h])) S->log_case = true;

    for (size_t h = 0; h < m; h++) {
        S->er[h] = Br[h]; S->ei[h] = Bi[h];
        /* Upper parameters 1 + b_h - a_k, lower 1 + b_h - b_j (j != h). */
        size_t nl = 0;
        for (size_t k = 0; k < p; k++) {
            S->ur[h][k] = 1.0 + Br[h] - Ar[k];
            S->ui[h][k] = Bi[h] - Ai[k];
        }
        for (size_t j = 0; j < q; j++) {
            if (j == h) continue;
            S->lr[h][nl] = 1.0 + Br[h] - Br[j];
            S->li[h][nl] = Bi[h] - Bi[j];
            nl++;
        }
        Prod C; prod_init(&C);
        /* Numerator Gammas Gamma(1 + b_h - a_j), j <= n. */
        for (size_t j = 0; j < n; j++)
            if (!prod_gamma(&C, 1.0 + Br[h] - Ar[j], Bi[h] - Ai[j], false)) return false;
        if (S->log_case) { S->cr[h] = 0.0; S->ci[h] = 0.0; continue; }
        /* 1/Gamma(a_j - b_h), j > n. */
        for (size_t j = n; j < p; j++)
            if (!prod_gamma(&C, Ar[j] - Br[h], Ai[j] - Bi[h], true)) return false;
        /* A lower parameter 1 + b_h - b_j (necessarily j > m: j <= m would be
         * the log case) that is a non-positive integer needs the regularized
         * form; the common case does not, and keeps the textbook coefficient. */
        bool shift = false;
        for (size_t l = 0; l < nl; l++)
            if (is_nonpos_int(S->lr[h][l], S->li[h][l])) shift = true;
        if (!shift) {
            /* Gamma(b_j - b_h), j <= m, j != h; 1/Gamma(1 + b_h - b_j), j > m:
             * the latter pairs with the term's own lower parameter. */
            for (size_t j = 0; j < q; j++) {
                if (j == h) continue;
                bool ok = (j < m) ? prod_gamma(&C, Br[j] - Br[h], Bi[j] - Bi[h], false)
                                  : prod_gamma(&C, 1.0 + Br[h] - Br[j], Bi[h] - Bi[j], true);
                if (!ok) return false;
            }
            S->kx[h] = 0;
            if (!prod_value(&C, &S->cr[h], &S->ci[h])) return false;
            if (S->ci[h] != 0.0) S->real = false;
            continue;
        }
        /* Regularized route. pi / sin(pi (b_j - b_h)), j <= m, j != h, is the
         * Gamma pair Gamma(b_j - b_h) Gamma(1 + b_h - b_j); the term is then
         * C * R(u; l; x) with R the REGULARIZED pFq, every lower parameter's
         * 1/Gamma applied below. */
        for (size_t j = 0; j < m; j++) {
            if (j == h) continue;
            double sr, si;
            c_sinpi(Br[j] - Br[h], Bi[j] - Bi[h], &sr, &si);
            prod_mul(&C, M_PI, 0.0);
            prod_div(&C, sr, si);
        }
        /* A non-positive integer lower parameter -N is shifted out exactly:
         *   R(u; -N, l'; x) = prod_k (u_k)_{N+1} x^{N+1} R(u+N+1; N+2, l'+N+1; x)
         * repeated until no such parameter remains. */
        int kx = 0;
        for (int guard = 0; guard < 64; guard++) {
            size_t hit = (size_t)-1;
            for (size_t l = 0; l < nl; l++)
                if (is_nonpos_int(S->lr[h][l], S->li[h][l])) { hit = l; break; }
            if (hit == (size_t)-1) break;
            /* bound before converting: out-of-range -> long is UB (G10) */
            if (!(-S->lr[h][hit] <= 1.0e6)) return false;
            long N = (long)(-S->lr[h][hit]);
            for (size_t k = 0; k < p; k++) {
                for (long t = 0; t <= N; t++) prod_mul(&C, S->ur[h][k] + (double)t, S->ui[h][k]);
                S->ur[h][k] += (double)(N + 1);
            }
            for (size_t l = 0; l < nl; l++) {
                if (l == hit) S->lr[h][l] = (double)(N + 2);
                else          S->lr[h][l] += (double)(N + 1);
            }
            kx += (int)(N + 1);
        }
        S->kx[h] = kx;
        /* 1/Gamma of every (now pole-free) lower parameter. */
        for (size_t l = 0; l < nl; l++)
            if (!prod_gamma(&C, S->lr[h][l], S->li[h][l], true)) return false;
        if (!prod_value(&C, &S->cr[h], &S->ci[h])) return false;
        if (S->ci[h] != 0.0) S->real = false;
    }
    S->valid = true;
    return true;
}

/* ------------------------------------------------------------------ */
/*  The log case: the residue table                                    */
/* ------------------------------------------------------------------ */

/* One linear factor of R(beta + k + eps): (c0 + dir k) + sg eps, in the
 * numerator (num) or the denominator. iz: c0 is an exact (snapped) integer, so
 * the constant term is zero exactly when c0 + dir k == 0. */
typedef struct { scdd c0; bool iz, num; double dir, sg; } MgLinFac;

/* A parameter as an exact complex double-double: the double (hi) plus the low
 * part the dual's 1 - x lost to rounding (lo; NULL for the primary
 * expansion, whose parameters are exactly the doubles). */
typedef struct { const double *r, *i, *rl, *il; } MgParams;

static scdd mg_param(const MgParams* X, size_t j) {
    scdd v;
    v.re.hi = X->r[j]; v.re.lo = X->rl ? X->rl[j] : 0.0;
    v.im.hi = X->i[j]; v.im.lo = X->il ? X->il[j] : 0.0;
    return v;
}

/* param - beta: the exact integer when the double difference is one (the rule
 * that detected the log case), else the exact double-double difference. The
 * low parts matter: the G-function can amplify a parameter's rounding by ~10^2
 * (MeijerG[{{1,2},{}},{{1/3,1/5},{}},5] moves 7e-15 under the 5.5e-17 that
 * 1 - b loses in double), so the dual's parameters are carried exactly. */
static scdd mg_diff_snap(scdd x, scdd beta, bool* iz) {
    double dr = x.re.hi - beta.re.hi, di = x.im.hi - beta.im.hi;
    if (is_int(dr, di)) { *iz = true; return scdd_from(dr, 0.0); }
    *iz = false;
    scdd d; d.re = sdd_sub(x.re, beta.re); d.im = sdd_sub(x.im, beta.im);
    return d;
}

static scdd mg_one_minus(scdd d) {
    scdd r; r.re = sdd_sub(sdd_from(1.0), d.re); r.im = sdd_neg(d.im); return r;
}

static scdd mg_minus_one(scdd d) { d.re = sdd_add_d(d.re, -1.0); return d; }

static void mg_put_x(double* x, scdd d) {
    x[0] = d.re.hi; x[1] = d.re.lo; x[2] = d.im.hi; x[3] = d.im.lo;
}

static void log_tab_free(SfMgLogTab* t) {
    if (!t) return;
    for (size_t c = 0; c < t->nc; c++) free(t->cl[c].c);
    free(t);
}

/* Tabulate one class with base beta and T members (see the header). */
static bool build_log_class(SfMgLogClass* C, scdd beta, size_t T,
                            const MgParams* A, size_t n, size_t p,
                            const MgParams* B, size_t m, size_t q) {
    double x[4 * 2 * SF_MG_MAXP];
    int sg[2 * SF_MG_MAXP], ex[2 * SF_MG_MAXP];
    MgLinFac lf[2 * SF_MG_MAXP];
    size_t nf = 0;
    /* F(beta + eps) = prod_{j<m} Gamma(D_j - eps) prod_{j<n} Gamma(1 - A_j + eps)
     *               / prod_{j>=m} Gamma(1 - D_j + eps) / prod_{j>=n} Gamma(A_j - eps),
     * with D_j = b_j - beta, A_j = a_j - beta. */
    for (size_t j = 0; j < q; j++, nf++) {
        bool iz;
        scdd D = mg_diff_snap(mg_param(B, j), beta, &iz);
        if (j < m) {
            mg_put_x(x + 4 * nf, D); sg[nf] = -1; ex[nf] = +1;
            lf[nf].c0 = mg_minus_one(D); lf[nf].dir = -1.0; lf[nf].sg = -1.0;
        } else {
            mg_put_x(x + 4 * nf, mg_one_minus(D)); sg[nf] = +1; ex[nf] = -1;
            lf[nf].c0 = mg_one_minus(D); lf[nf].dir = +1.0; lf[nf].sg = +1.0;
        }
        lf[nf].iz = iz; lf[nf].num = false;
    }
    for (size_t j = 0; j < p; j++, nf++) {
        bool iz;
        scdd D = mg_diff_snap(mg_param(A, j), beta, &iz);
        if (j < n) {
            /* Gamma(1 - A_j + s) would have a pole on the class's poles: the
             * contour cannot separate them (an inconsistent parameter set). */
            if (iz && D.re.hi >= 1.0) return false;
            mg_put_x(x + 4 * nf, mg_one_minus(D)); sg[nf] = +1; ex[nf] = +1;
            lf[nf].c0 = mg_one_minus(D); lf[nf].dir = +1.0; lf[nf].sg = +1.0;
        } else {
            mg_put_x(x + 4 * nf, D); sg[nf] = -1; ex[nf] = -1;
            lf[nf].c0 = mg_minus_one(D); lf[nf].dir = -1.0; lf[nf].sg = -1.0;
        }
        lf[nf].iz = iz; lf[nf].num = true;
    }

    double f0[4 * SF_MG_MAXP];
    long v;
    if (!flint_mg_laurent(nf, x, sg, ex, T, &v, f0)) return false;
    scdd f[SF_MG_MAXP];
    for (size_t i = 0; i < T; i++) {
        f[i].re.hi = f0[4 * i];     f[i].re.lo = f0[4 * i + 1];
        f[i].im.hi = f0[4 * i + 2]; f[i].im.lo = f0[4 * i + 3];
    }
    C->br = beta.re.hi + beta.re.lo; C->bi = beta.im.hi + beta.im.lo; C->T = T;
    C->c = calloc((size_t)SF_MG_LOG_KMAX * T * 4, sizeof(double));
    if (!C->c) return false;

    double fact = 1.0;                       /* j!, exact for j < 18 */
    double ifact[SF_MG_MAXP];
    for (size_t j = 0; j < T; j++) { if (j) fact *= (double)j; ifact[j] = fact; }

    size_t K = 0;
    for (size_t k = 0; k < SF_MG_LOG_KMAX; k++) {
        if (-v > (long)T) return false;      /* pole order beyond the class: inconsistent */
        /* Residue row: c[k][j] = -f_{-1-j} / j!, f_{e} at index e - v. */
        double* row = C->c + 4 * k * T;
        bool tiny = true;
        for (size_t j = 0; j < T; j++) {
            long idx = -1 - (long)j - v;
            if (idx < 0 || idx >= (long)T) continue;
            scdd cc;
            cc.re = sdd_neg(sdd_div_d(f[idx].re, ifact[j]));
            cc.im = sdd_neg(sdd_div_d(f[idx].im, ifact[j]));
            row[4 * j] = cc.re.hi; row[4 * j + 1] = cc.re.lo;
            row[4 * j + 2] = cc.im.hi; row[4 * j + 3] = cc.im.lo;
        }
        K = k + 1;
        for (size_t i = 0; i < T; i++) {
            if (!isfinite(f[i].re.hi) || !isfinite(f[i].im.hi)) return false;
            if (fabs(f[i].re.hi) > 1e-280 || fabs(f[i].im.hi) > 1e-280) tiny = false;
        }
        if (tiny) break;                     /* every later row underflows too */
        /* f <- f * R(beta + k + eps). */
        for (size_t i = 0; i < nf; i++) {
            scdd c = lf[i].c0;
            c.re = sdd_add_d(c.re, lf[i].dir * (double)k);
            if (lf[i].iz && c.re.hi == 0.0) {
                /* (sg eps)^{+-1}: the pole order moves, the coefficients only flip sign. */
                v += lf[i].num ? 1 : -1;
                if (lf[i].sg < 0.0)
                    for (size_t t = 0; t < T; t++) { f[t].re = sdd_neg(f[t].re); f[t].im = sdd_neg(f[t].im); }
                continue;
            }
            if (lf[i].num) {
                for (size_t t = T; t-- > 1;)
                    f[t] = scdd_add(scdd_mul(c, f[t]), scdd_mul_d(f[t - 1], lf[i].sg));
                f[0] = scdd_mul(c, f[0]);
            } else {
                scdd ic = scdd_inv(c);
                f[0] = scdd_mul(f[0], ic);
                for (size_t t = 1; t < T; t++)
                    f[t] = scdd_mul(scdd_add(f[t], scdd_mul_d(f[t - 1], -lf[i].sg)), ic);
            }
        }
    }
    C->K = K;
    return true;
}

/* Partition b_1..b_m into integer-spaced classes and tabulate each. NULL when
 * a class cannot be tabulated (no FLINT, an inconsistent set, overflow): the
 * series then declines, as it did before the table existed. */
static SfMgLogTab* build_log_table(const MgParams* A, size_t n, size_t p,
                                   const MgParams* B, size_t m, size_t q) {
    const double *Br = B->r, *Bi = B->i;
    SfMgLogTab* t = calloc(1, sizeof *t);
    if (!t) return NULL;
    int cls[SF_MG_MAXP];
    for (size_t j = 0; j < m; j++) cls[j] = -1;
    for (size_t h = 0; h < m; h++) {
        if (cls[h] >= 0) continue;
        int c = (int)t->nc;
        size_t T = 0, base = h;
        double lo = 0.0;
        for (size_t j = h; j < m; j++) {
            if (cls[j] >= 0 || !is_int(Br[j] - Br[h], Bi[j] - Bi[h])) continue;
            cls[j] = c; T++;
            if (Br[j] - Br[h] < lo) { lo = Br[j] - Br[h]; base = j; }
        }
        scdd beta = mg_param(B, base);
        double br = Br[base], bi = Bi[base];
        for (size_t j = h; j < m; j++)        /* integer offsets >= 0 from the base */
            if (cls[j] == c && (!is_int(Br[j] - br, Bi[j] - bi) || Br[j] - br < 0.0)) {
                log_tab_free(t); return NULL;
            }
        t->nc++;
        if (!build_log_class(&t->cl[c], beta, T, A, n, p, B, m, q)) {
            log_tab_free(t); return NULL;
        }
    }
    return t;
}

void sf_mg_plan_free(SfMgPlan* P) {
    if (!P) return;
    log_tab_free(P->zs.lt); P->zs.lt = NULL;
    log_tab_free(P->ds.lt); P->ds.lt = NULL;
}

bool sf_mg_plan(const double* ar, const double* ai, size_t n, size_t p,
                const double* br, const double* bi, size_t m, size_t q,
                double r, SfMgPlan* P) {
    memset(P, 0, sizeof *P);
    if (p > SF_MG_MAXP || q > SF_MG_MAXP || n > p || m > q) return false;
    if (!isfinite(r) || r == 0.0) return false;
    P->n = n; P->p = p; P->m = m; P->q = q; P->r = r;
    P->real_params = true;
    for (size_t j = 0; j < p; j++) {
        P->ar[j] = ar[j]; P->ai[j] = ai ? ai[j] : 0.0;
        if (!isfinite(P->ar[j]) || !isfinite(P->ai[j])) return false;
        if (P->ai[j] != 0.0) P->real_params = false;
    }
    for (size_t j = 0; j < q; j++) {
        P->br[j] = br[j]; P->bi[j] = bi ? bi[j] : 0.0;
        if (!isfinite(P->br[j]) || !isfinite(P->bi[j])) return false;
        if (P->bi[j] != 0.0) P->real_params = false;
    }
    if (p <= q) {
        if (!build_series(&P->zs, P->ar, P->ai, n, p, P->br, P->bi, m, q))
            P->zs.valid = false;
        else if (P->zs.log_case) {
            MgParams A = { P->ar, P->ai, NULL, NULL }, B = { P->br, P->bi, NULL, NULL };
            P->zs.lt = build_log_table(&A, n, p, &B, m, q);
        }
    }
    if (p >= q) {
        /* Dual: G^{n,m}_{q,p}(1 - b; 1 - a; 1/w). */
        double Ar[SF_MG_MAXP], Ai[SF_MG_MAXP], Br[SF_MG_MAXP], Bi[SF_MG_MAXP];
        for (size_t j = 0; j < q; j++) { Ar[j] = 1.0 - P->br[j]; Ai[j] = -P->bi[j]; }
        for (size_t j = 0; j < p; j++) { Br[j] = 1.0 - P->ar[j]; Bi[j] = -P->ai[j]; }
        if (!build_series(&P->ds, Ar, Ai, m, q, Br, Bi, n, p)) P->ds.valid = false;
        else if (P->ds.log_case) {
            /* 1 - x exactly: the rounded double plus what it lost. */
            double Arl[SF_MG_MAXP], Brl[SF_MG_MAXP];
            for (size_t j = 0; j < q; j++) Arl[j] = sdd_two_sum(1.0, -P->br[j]).lo;
            for (size_t j = 0; j < p; j++) Brl[j] = sdd_two_sum(1.0, -P->ar[j]).lo;
            MgParams A = { Ar, Ai, Arl, NULL }, B = { Br, Bi, Brl, NULL };
            P->ds.lt = build_log_table(&A, m, q, &B, n, p);
        }
    }
    P->ok = true;
    return true;
}

/* ------------------------------------------------------------------ */
/*  Evaluation                                                         */
/* ------------------------------------------------------------------ */

/* pFq(u; l; x) for one term, real lane. */
static bool term_pfq_r(const SfMgSeries* S, size_t h, double x, double* F) {
    if (S->pu == 0 && S->ql == 0) { *F = exp(x); return isfinite(*F); }
    return sf_machine_pfq(S->ur[h], S->pu, S->lr[h], S->ql, x, F);
}

/* pFq(u; l; x) for one term, complex lane. */
static bool term_pfq_c(const SfMgSeries* S, size_t h, double xr, double xi,
                       double* Fr, double* Fi) {
    if (S->pu == 0 && S->ql == 0) { c_exp(xr, xi, Fr, Fi); return isfinite(*Fr) && isfinite(*Fi); }
    if (S->pu == 1 && S->ql == 0) {
        /* 1F0(u;;x) = (1 - x)^(-u), principal branch: the closed form is valid
         * on the whole plane bar the cut, so |x| = 1 needs no series. */
        if (xr == 1.0 && xi == 0.0) return false;
        double lr, li, er, ei;
        c_log(1.0 - xr, -xi, &lr, &li);
        c_mul(-S->ur[h][0], -S->ui[h][0], lr, li, &er, &ei);
        c_exp(er, ei, Fr, Fi);
        return isfinite(*Fr) && isfinite(*Fi);
    }
    return sf_machine_pfq_c(S->ur[h], S->ui[h], S->pu, S->lr[h], S->li[h], S->ql,
                            xr, xi, Fr, Fi);
}

/* Sum one expansion at log w = (lr, li). */
static bool eval_series(const SfMgSeries* S, double lr, double li, bool real_lane,
                        double* rr, double* ri) {
    if (!S->valid || S->log_case) return false;
    if (S->nt == 0) { *rr = 0.0; *ri = 0.0; return true; }
    double wr, wi;
    c_exp(lr, li, &wr, &wi);
    double xr = S->sgn * wr, xi = S->sgn * wi;
    double sr = 0.0, si = 0.0, maxabs = 0.0;
    for (size_t h = 0; h < S->nt; h++) {
        if (S->cr[h] == 0.0 && S->ci[h] == 0.0) continue;
        double tr, ti;
        if (real_lane) {
            double F;
            if (!term_pfq_r(S, h, xr, &F)) return false;
            /* w^{e} for real w > 0 and real e. */
            double pw = exp(S->er[h] * lr);
            tr = S->cr[h] * pw * F;
            for (int k = 0; k < S->kx[h]; k++) tr *= xr;
            ti = 0.0;
        } else {
            double Fr, Fi, er, ei, pr, pi;
            if (!term_pfq_c(S, h, xr, xi, &Fr, &Fi)) return false;
            c_mul(S->er[h], S->ei[h], lr, li, &er, &ei);
            c_exp(er, ei, &pr, &pi);
            c_mul(S->cr[h], S->ci[h], pr, pi, &tr, &ti);
            c_mul(tr, ti, Fr, Fi, &tr, &ti);
            for (int k = 0; k < S->kx[h]; k++) c_mul(tr, ti, xr, xi, &tr, &ti);
        }
        if (!isfinite(tr) || !isfinite(ti)) return false;
        double a = c_abs(tr, ti);
        if (a > maxabs) maxabs = a;
        sr += tr; si += ti;
    }
    double s = c_abs(sr, si);
    if (maxabs > 0.0) {
        if (s == 0.0) return false;                     /* total cancellation */
        if (log2(maxabs / s) > SF_MG_MAX_LOST_BITS) return false;
    }
    *rr = sr; *ri = si;
    return true;
}

/* ---- the log case, per point ---------------------------------------------- */

/* Termination: stop after three consecutive terms below tol times the largest
 * (the series is entire for p < q, geometric for p == q, so once past the peak
 * the terms only shrink); 2^-60 for the double sum, 2^-112 for the dd one. */
#define MG_LOG_TOL_D  8.673617379884035e-19
#define MG_LOG_TOL_DD 1.925929944387236e-34

/* Each class sum V = sum_j (sum_k c[k][j] w^k) L^j, with its term mass
 * M = sum_k |w|^k sum_j |c[k][j]| |L|^j (the cancellation measure). Returns
 * false when the table ends before the terms do. */
static bool lc_real_d(const SfMgLogClass* C, double w, double L, double* V, double* M) {
    size_t T = C->T;
    double S[SF_MG_MAXP], lp[SF_MG_MAXP];
    for (size_t j = 0; j < T; j++) { S[j] = 0.0; lp[j] = j ? lp[j - 1] * fabs(L) : 1.0; }
    double wk = 1.0, mass = 0.0, amax = 0.0;
    int small = 0;
    bool conv = false;
    const double* c = C->c;
    for (size_t k = 0; k < C->K; k++, c += 4 * T) {
        double a = 0.0;
        for (size_t j = 0; j < T; j++) { S[j] += c[4 * j] * wk; a += fabs(c[4 * j]) * lp[j]; }
        a *= fabs(wk);
        mass += a;
        if (a > amax) amax = a;
        if (amax > 0.0 && a <= MG_LOG_TOL_D * amax) { if (++small >= 3) { conv = true; break; } }
        else small = 0;
        wk *= w;
    }
    if (!conv && amax > 0.0) return false;
    double v = S[T - 1];
    for (size_t j = T - 1; j-- > 0;) v = v * L + S[j];
    *V = v; *M = mass;
    return isfinite(v) && isfinite(mass);
}

static bool lc_cplx_d(const SfMgLogClass* C, double wr, double wi, double Lr, double Li,
                      double* Vr, double* Vi, double* M) {
    size_t T = C->T;
    double Sr[SF_MG_MAXP], Si[SF_MG_MAXP], lp[SF_MG_MAXP];
    double aL = c_abs(Lr, Li), aw = c_abs(wr, wi);
    for (size_t j = 0; j < T; j++) { Sr[j] = Si[j] = 0.0; lp[j] = j ? lp[j - 1] * aL : 1.0; }
    double pr = 1.0, pi = 0.0, pm = 1.0, mass = 0.0, amax = 0.0;
    int small = 0;
    bool conv = false;
    const double* c = C->c;
    for (size_t k = 0; k < C->K; k++, c += 4 * T) {
        double a = 0.0;
        for (size_t j = 0; j < T; j++) {
            double cr = c[4 * j], ci = c[4 * j + 2];
            Sr[j] += cr * pr - ci * pi;
            Si[j] += cr * pi + ci * pr;
            a += (fabs(cr) + fabs(ci)) * lp[j];
        }
        a *= pm;
        mass += a;
        if (a > amax) amax = a;
        if (amax > 0.0 && a <= MG_LOG_TOL_D * amax) { if (++small >= 3) { conv = true; break; } }
        else small = 0;
        c_mul(pr, pi, wr, wi, &pr, &pi);
        pm *= aw;
    }
    if (!conv && amax > 0.0) return false;
    double vr = Sr[T - 1], vi = Si[T - 1];
    for (size_t j = T - 1; j-- > 0;) {
        c_mul(vr, vi, Lr, Li, &vr, &vi);
        vr += Sr[j]; vi += Si[j];
    }
    *Vr = vr; *Vi = vi; *M = mass;
    return isfinite(vr) && isfinite(vi) && isfinite(mass);
}

static bool lc_real_dd(const SfMgLogClass* C, sdd w, sdd L, double* V, double* M) {
    size_t T = C->T;
    sdd S[SF_MG_MAXP];
    double lp[SF_MG_MAXP];
    for (size_t j = 0; j < T; j++) { S[j] = sdd_from(0.0); lp[j] = j ? lp[j - 1] * fabs(L.hi) : 1.0; }
    sdd wk = sdd_from(1.0);
    double mass = 0.0, amax = 0.0;
    int small = 0;
    bool conv = false;
    const double* c = C->c;
    for (size_t k = 0; k < C->K; k++, c += 4 * T) {
        double a = 0.0;
        for (size_t j = 0; j < T; j++) {
            sdd cj = { c[4 * j], c[4 * j + 1] };
            S[j] = sdd_add(S[j], sdd_mul(cj, wk));
            a += fabs(cj.hi) * lp[j];
        }
        a *= fabs(wk.hi);
        mass += a;
        if (a > amax) amax = a;
        if (amax > 0.0 && a <= MG_LOG_TOL_DD * amax) { if (++small >= 3) { conv = true; break; } }
        else small = 0;
        wk = sdd_mul(wk, w);
    }
    if (!conv && amax > 0.0) return false;
    sdd v = S[T - 1];
    for (size_t j = T - 1; j-- > 0;) v = sdd_add(sdd_mul(v, L), S[j]);
    *V = v.hi + v.lo; *M = mass;
    return isfinite(*V) && isfinite(mass);
}

static bool lc_cplx_dd(const SfMgLogClass* C, scdd w, scdd L, double* Vr, double* Vi, double* M) {
    size_t T = C->T;
    scdd S[SF_MG_MAXP];
    double lp[SF_MG_MAXP];
    double aL = c_abs(L.re.hi, L.im.hi), aw = c_abs(w.re.hi, w.im.hi);
    for (size_t j = 0; j < T; j++) { S[j] = scdd_from(0.0, 0.0); lp[j] = j ? lp[j - 1] * aL : 1.0; }
    scdd wk = scdd_from(1.0, 0.0);
    double pm = 1.0, mass = 0.0, amax = 0.0;
    int small = 0;
    bool conv = false;
    const double* c = C->c;
    for (size_t k = 0; k < C->K; k++, c += 4 * T) {
        double a = 0.0;
        for (size_t j = 0; j < T; j++) {
            scdd cj; cj.re.hi = c[4 * j]; cj.re.lo = c[4 * j + 1];
            cj.im.hi = c[4 * j + 2]; cj.im.lo = c[4 * j + 3];
            S[j] = scdd_add(S[j], scdd_mul(cj, wk));
            a += (fabs(cj.re.hi) + fabs(cj.im.hi)) * lp[j];
        }
        a *= pm;
        mass += a;
        if (a > amax) amax = a;
        if (amax > 0.0 && a <= MG_LOG_TOL_DD * amax) { if (++small >= 3) { conv = true; break; } }
        else small = 0;
        wk = scdd_mul(wk, w);
        pm *= aw;
    }
    if (!conv && amax > 0.0) return false;
    scdd v = S[T - 1];
    for (size_t j = T - 1; j-- > 0;) v = scdd_add(scdd_mul(v, L), S[j]);
    *Vr = v.re.hi + v.re.lo; *Vi = v.im.hi + v.im.lo; *M = mass;
    return isfinite(*Vr) && isfinite(*Vi) && isfinite(mass);
}

/* The log-case series at the point z (zr, zi) of a plan with root r, sl = +1
 * for the expansion in w, -1 for the one in 1/w; (lr, li) = sl log(z) / r in
 * double. Each class is summed in double first and re-summed in double-double
 * only when its measured cancellation exceeds SF_MG_LOG_D_LOST_BITS (the double
 * sum's error is ~ 2^lost ulp, so this keeps it to a few ulp). */
static bool eval_log_series(const SfMgSeries* S, double zr, double zi, double r, double sl,
                            double lr, double li, bool real_lane, double* rr, double* ri) {
    const SfMgLogTab* t = S->lt;
    if (!S->valid || !t) return false;
    const double lim_d = ldexp(1.0, (int)SF_MG_LOG_D_LOST_BITS);
    const double lim_dd = ldexp(1.0, (int)SF_MG_DD_MAX_LOST_BITS);
    bool direct = (r == 1.0 && sl > 0.0);              /* w = z exactly */
    double wr, wi;
    if (direct) { wr = zr; wi = zi; } else c_exp(lr, li, &wr, &wi);
    if (real_lane) wi = 0.0;
    bool have_dd = false;
    sdd wd = sdd_from(0.0), Ld = sdd_from(0.0);
    scdd wc = scdd_from(0.0, 0.0), Lc = scdd_from(0.0, 0.0);
    double sr = 0.0, si = 0.0, smass = 0.0;
    for (size_t c = 0; c < t->nc; c++) {
        const SfMgLogClass* C = &t->cl[c];
        double Vr = 0.0, Vi = 0.0, M = 0.0;
        bool ok = real_lane ? lc_real_d(C, wr, lr, &Vr, &M)
                            : lc_cplx_d(C, wr, wi, lr, li, &Vr, &Vi, &M);
        if (!(ok && (M == 0.0 || M <= c_abs(Vr, Vi) * lim_d))) {
            if (!have_dd) {
                if (real_lane) {
                    Ld = sdd_log(sdd_from(zr));
                    if (r != 1.0) Ld = sdd_div_d(Ld, r);
                    if (sl < 0.0) Ld = sdd_neg(Ld);
                    wd = direct ? sdd_from(zr)
                       : (r == 1.0 ? sdd_div(sdd_from(1.0), sdd_from(zr)) : sdd_exp(Ld));
                } else {
                    Lc = scdd_log_d(zr, zi);
                    if (r != 1.0) { Lc.re = sdd_div_d(Lc.re, r); Lc.im = sdd_div_d(Lc.im, r); }
                    if (sl < 0.0) { Lc.re = sdd_neg(Lc.re); Lc.im = sdd_neg(Lc.im); }
                    wc = direct ? scdd_from(zr, zi) : scdd_exp(Lc);
                }
                have_dd = true;
            }
            ok = real_lane ? lc_real_dd(C, wd, Ld, &Vr, &M)
                           : lc_cplx_dd(C, wc, Lc, &Vr, &Vi, &M);
            if (!ok || (M > 0.0 && !(M <= c_abs(Vr, Vi) * lim_dd))) return false;
        }
        /* times w^beta = exp(beta log w) */
        double er, ei, pr, pi, tr, ti;
        c_mul(C->br, C->bi, lr, li, &er, &ei);
        c_exp(er, ei, &pr, &pi);
        c_mul(Vr, Vi, pr, pi, &tr, &ti);
        if (!isfinite(tr) || !isfinite(ti)) return false;
        sr += tr; si += ti;
        smass += c_abs(tr, ti);
    }
    double s = c_abs(sr, si);
    if (smass > 0.0) {
        if (s == 0.0) return false;                     /* total cancellation */
        if (log2(smass / s) > SF_MG_MAX_LOST_BITS) return false;
    }
    *rr = sr; *ri = real_lane ? 0.0 : si;
    return true;
}

bool sf_mg_eval_machine(const SfMgPlan* P, double zr, double zi,
                        double* rr, double* ri, bool* is_c) {
    if (!P->ok || !isfinite(zr) || !isfinite(zi)) return false;
    if (zr == 0.0 && zi == 0.0) return false;           /* the exact path owns z = 0 */
    double lr, li;
    c_log(zr, zi, &lr, &li);
    lr /= P->r; li /= P->r;
    const SfMgSeries* S;
    double sl = 1.0;                                    /* +1: w, -1: 1/w */
    if (P->p < P->q) S = &P->zs;
    else if (P->p > P->q) { S = &P->ds; sl = -1.0; }
    else {
        if (lr < 0.0) S = &P->zs;
        else if (lr > 0.0) { S = &P->ds; sl = -1.0; }
        else {
            /* |w| = 1 exactly: only the shapes whose term function continues
             * analytically past the unit circle (1F0, 2F1). */
            if (P->p > 2) return false;
            S = &P->zs;
        }
    }
    bool real_lane = P->real_params && S->real && li == 0.0 && zi == 0.0 && zr > 0.0;
    if (S->log_case) {
        if (!eval_log_series(S, zr, zi, P->r, sl, sl * lr, real_lane ? 0.0 : sl * li,
                             real_lane, rr, ri))
            return false;
        *is_c = !real_lane;
        return true;
    }
    if (real_lane && eval_series(S, sl * lr, 0.0, true, rr, ri)) {
        *ri = 0.0; *is_c = false; return true;
    }
    if (!eval_series(S, sl * lr, sl * li, false, rr, ri)) return false;
    *is_c = true;
    if (real_lane) { *ri = 0.0; *is_c = false; }       /* real by construction */
    return true;
}

bool sf_mg_eval(const SfMgPlan* P, double zr, double zi,
                double* rr, double* ri, bool* is_c) {
    if (sf_mg_eval_machine(P, zr, zi, rr, ri, is_c)) return true;
    if (!P->ok || !isfinite(zr) || !isfinite(zi)) return false;
    if (zr == 0.0 && zi == 0.0) return false;
    return flint_mg_eval_d(P->ar, P->ai, P->n, P->p, P->br, P->bi, P->m, P->q,
                           P->r, zr, zi, rr, ri, is_c);
}

/* ------------------------------------------------------------------ */
/*  Flattened n-ary layout (ndkernels.c, Compile[])                    */
/* ------------------------------------------------------------------ */

double sf_mg_shape_code(size_t n, size_t p, size_t m) {
    return (double)n + 64.0 * (double)p + 4096.0 * (double)m;
}

/* A small, append-only cache of plans keyed by the parameter block, so a
 * compiled loop over z builds the plan once. Entries are never replaced or
 * freed while the process runs (plans are immutable once published), which is
 * what makes the lookup safe from the threaded fused-map loops without holding
 * the lock during evaluation; once the cache is full, further parameter sets
 * build their plan per call on the stack. */
#define MG_MEMO_SLOTS 16
#define MG_KEY_MAX (2 + 4 * SF_MG_MAXP)
typedef struct { size_t len; double key[MG_KEY_MAX]; SfMgPlan* plan; } MgMemo;
static MgMemo g_memo[MG_MEMO_SLOTS];
static size_t g_memo_n = 0;
#ifdef MATHILDA_THREADS
static pthread_mutex_t g_memo_lock = PTHREAD_MUTEX_INITIALIZER;
#define MG_LOCK()   pthread_mutex_lock(&g_memo_lock)
#define MG_UNLOCK() pthread_mutex_unlock(&g_memo_lock)
#else
#define MG_LOCK()   ((void)0)
#define MG_UNLOCK() ((void)0)
#endif

static const SfMgPlan* memo_find(const double* key, size_t len) {
    const SfMgPlan* hit = NULL;
    MG_LOCK();
    for (size_t s = 0; s < g_memo_n; s++)
        if (g_memo[s].len == len && memcmp(g_memo[s].key, key, len * sizeof(double)) == 0) {
            hit = g_memo[s].plan; break;
        }
    MG_UNLOCK();
    return hit;
}

/* Publish P (moving its heap tables into the memo's copy) and return the
 * published plan, or NULL when the memo is full (P keeps its tables). */
static const SfMgPlan* memo_add(const double* key, size_t len, SfMgPlan* P) {
    SfMgPlan* c = NULL;
    MG_LOCK();
    if (g_memo_n < MG_MEMO_SLOTS) {
        c = malloc(sizeof *c);
        if (c) {
            *c = *P;
            P->zs.lt = NULL; P->ds.lt = NULL;        /* ownership moved */
            g_memo[g_memo_n].len = len;
            memcpy(g_memo[g_memo_n].key, key, len * sizeof(double));
            g_memo[g_memo_n].plan = c;
            g_memo_n++;
        }
    }
    MG_UNLOCK();
    return c;
}

const SfMgPlan* sf_mg_plan_get(const double* ar, const double* ai, size_t n, size_t p,
                               const double* br, const double* bi, size_t m, size_t q,
                               double r, bool memo_all, SfMgPlan* scratch) {
    memset(scratch, 0, sizeof *scratch);
    if (p > SF_MG_MAXP || q > SF_MG_MAXP || n > p || m > q) return NULL;
    /* Key: shape code, r, then (re, im) of every parameter; q is implied by
     * the length. */
    double key[MG_KEY_MAX];
    size_t len = 0;
    key[len++] = sf_mg_shape_code(n, p, m); key[len++] = r;
    for (size_t j = 0; j < p; j++) { key[len++] = ar[j]; key[len++] = ai ? ai[j] : 0.0; }
    for (size_t j = 0; j < q; j++) { key[len++] = br[j]; key[len++] = bi ? bi[j] : 0.0; }
    const SfMgPlan* P = memo_find(key, len);
    if (P) return P;
    if (!sf_mg_plan(ar, ai, n, p, br, bi, m, q, r, scratch)) {
        sf_mg_plan_free(scratch);
        return NULL;
    }
    if (memo_all || scratch->zs.lt || scratch->ds.lt) {
        const SfMgPlan* pub = memo_add(key, len, scratch);
        if (pub) return pub;
    }
    return scratch;
}

bool sf_mg_eval_flat(const double* re, const double* im, size_t count,
                     double* rr, double* ri) {
    if (count < 3 || im[0] != 0.0 || im[1] != 0.0) return false;
    double code = re[0];
    if (!(code >= 0.0) || code != nearbyint(code)) return false;
    long c = (long)code;
    size_t n = (size_t)(c % 64), p = (size_t)((c / 64) % 64), m = (size_t)(c / 4096);
    if (count < p + 3) return false;
    size_t q = count - p - 3;
    if (p > SF_MG_MAXP || q > SF_MG_MAXP) return false;
    SfMgPlan local;
    const SfMgPlan* P = sf_mg_plan_get(re + 2, im + 2, n, p, re + 2 + p, im + 2 + p, m, q,
                                       re[1], true, &local);
    if (!P) return false;
    bool is_c;
    bool ok = sf_mg_eval(P, re[count - 1], im[count - 1], rr, ri, &is_c);
    if (P == &local) sf_mg_plan_free(&local);
    return ok;
}

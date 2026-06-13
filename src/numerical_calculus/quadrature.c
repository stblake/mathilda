/*
 * quadrature.c — periodic-trapezoidal contour integration core.
 *
 * See quadrature.h for the public contract and the mathematical background.
 *
 * Layout:
 *   - shared helpers (penalty metric, decay classification);
 *   - machine path (double _Complex): trapezoidal sum, N-doubling with
 *     sample reuse + Aitken extrapolation, adaptive-radius search;
 *   - MPFR path: the same algorithm with a file-local complex toolkit
 *     (pairs of mpfr_t — no MPC library is linked).
 */

#include "quadrature.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Initial sample count for the machine path. Geometric convergence makes a
 * small N sufficient for a well-separated pole; doubling handles the rest. */
#define QD_N0_MACHINE 16

/* Decay ratio above which the doubling errors are deemed non-geometric
 * (i.e. the contour is probably crossing a branch cut / non-analytic). */
#define QD_DECAY_BRANCHCUT 0.7

/* Max relative jump between adjacent samples above which the integrand is
 * deemed discontinuous on the contour (a branch-cut crossing). At the final,
 * largest N a smooth integrand has tiny relative neighbour differences; a
 * sign-flipping branch cut keeps an O(1) jump no matter how fine the mesh. */
#define QD_REL_JUMP_BRANCHCUT 0.5

/* ------------------------------------------------------------------ *
 *  Shared: convergence classification + radius-search penalty         *
 * ------------------------------------------------------------------ */

/* Map a finished fixed-radius run to a scalar penalty (lower is better) for
 * the adaptive-radius search. The residue is radius-independent, so we just
 * reward converged runs that used few samples / reached a small error. */
static double qd_penalty(QdStatus status, double abs_err, int n_used) {
    switch (status) {
        case QD_OK:        return log10(abs_err + 1e-300) + 1e-4 * (double)n_used;
        case QD_NOCONV:    return 50.0 + 1e-4 * (double)n_used;
        case QD_BRANCHCUT: return 100.0;
        default:           return 1000.0;
    }
}

/* ================================================================== *
 *  Machine-precision path                                             *
 * ================================================================== */

/* Aitken / Shanks (Delta^2) extrapolation of three successive estimates of
 * a (complex) geometric sequence. Returns the accelerated value and sets
 * *err to a conservative gauge. Falls back to s2 when the second difference
 * is numerically negligible (already converged or noise-limited). */
static double _Complex qd_aitken_m(double _Complex s0, double _Complex s1,
                                   double _Complex s2, double* err) {
    double _Complex d1 = s1 - s0;
    double _Complex d2 = s2 - s1;
    double _Complex den = d2 - d1;
    double aden = cabs(den);
    if (aden <= 1e-300 || aden < 1e-12 * (cabs(d1) + cabs(d2))) {
        *err = cabs(d2);
        return s2;
    }
    double _Complex acc = s2 - (d2 * d2) / den;
    *err = cabs(acc - s2);
    return acc;
}

/* Trapezoidal residue estimate (r/N) sum f_k e^{i theta_k}. */
static double _Complex qd_trap_m(const double _Complex* f, int N, double r) {
    double _Complex acc = 0.0;
    for (int k = 0; k < N; k++) {
        double th = 2.0 * M_PI * (double)k / (double)N;
        acc += f[k] * (cos(th) + I * sin(th));
    }
    return acc * (r / (double)N);
}

/* Largest relative jump |f_{k+1}-f_k| / (|f_{k+1}|+|f_k|) around the ring.
 * A robust discontinuity gauge: smooth integrands shrink it as ~1/N, a
 * branch-cut sign flip keeps it O(1). */
static double qd_maxreljump_m(const double _Complex* f, int N) {
    double m = 0.0;
    for (int k = 0; k < N; k++) {
        double _Complex a = f[k], b = f[(k + 1) % N];
        double den = cabs(a) + cabs(b);
        if (den > 0.0 && isfinite(den)) {
            double rj = cabs(b - a) / den;
            if (rj > m) m = rj;
        }
    }
    return m;
}

/* Residue at a fixed radius via adaptive N-doubling with sample reuse. */
static QdResultMachine qd_fixed_m(QdSampleMachine f, void* ctx,
                                  double _Complex z0, double r,
                                  double prec_goal_digits, int max_recursion) {
    QdResultMachine R;
    R.value = 0.0; R.abs_err = INFINITY; R.status = QD_NONNUMERIC;
    R.n_used = 0; R.radius = r;

    double tol = pow(10.0, -prec_goal_digits);
    int N = QD_N0_MACHINE;

    double _Complex* fv = malloc(sizeof(double _Complex) * (size_t)N);
    if (!fv) return R;
    for (int k = 0; k < N; k++) {
        double th = 2.0 * M_PI * (double)k / (double)N;
        double _Complex z = z0 + r * (cos(th) + I * sin(th));
        if (!f(ctx, z, &fv[k])
            || !isfinite(creal(fv[k])) || !isfinite(cimag(fv[k]))) {
            free(fv); return R;   /* status stays QD_NONNUMERIC */
        }
    }

    double _Complex s_a = 0.0, s_b = 0.0, s_c = qd_trap_m(fv, N, r);
    int have = 1;
    double prev_err = INFINITY, decay_min = 1.0;
    double _Complex value = s_c;
    double abs_err = INFINITY;
    QdStatus st = QD_NOCONV;

    for (int it = 0; it < max_recursion; it++) {
        int newN = 2 * N;
        double _Complex* nf = realloc(fv, sizeof(double _Complex) * (size_t)newN);
        if (!nf) { free(fv); R.status = QD_NONNUMERIC; return R; }
        fv = nf;
        /* Old node k (angle 2pi k/N) maps to new even index 2k. Shift down
         * from the top so we never clobber a not-yet-moved entry. */
        for (int k = N - 1; k >= 0; k--) fv[2 * k] = fv[k];
        /* Sample the new odd nodes. */
        bool ok = true;
        for (int j = 0; j < N; j++) {
            int idx = 2 * j + 1;
            double th = 2.0 * M_PI * (double)idx / (double)newN;
            double _Complex z = z0 + r * (cos(th) + I * sin(th));
            if (!f(ctx, z, &fv[idx])
                || !isfinite(creal(fv[idx])) || !isfinite(cimag(fv[idx]))) {
                ok = false; break;
            }
        }
        if (!ok) { free(fv); R.status = QD_NONNUMERIC; return R; }
        N = newN;

        double _Complex S = qd_trap_m(fv, N, r);
        double err = cabs(S - s_c);
        double ratio = (isfinite(prev_err) && prev_err > 0.0) ? err / prev_err
                                                              : (err > 0.0 ? 1.0 : 0.0);
        if (ratio < decay_min) decay_min = ratio;
        prev_err = err;

        s_a = s_b; s_b = s_c; s_c = S; if (have < 3) have++;

        /* Best available estimate: the extrapolated value if it tightened
         * the error, else the raw doubling value. */
        double _Complex est = S;
        double est_err = err;
        if (have == 3) {
            double ae;
            double _Complex acc = qd_aitken_m(s_a, s_b, s_c, &ae);
            if (ae < est_err) { est = acc; est_err = ae; }
        }
        value = est; abs_err = est_err;

        /* Gate convergence on the RAW doubling error only — Aitken can make
         * `est_err` spuriously tiny for a non-geometric (discontinuous)
         * sequence, which would mask a branch cut. */
        double scale = 1.0 + cabs(s_c);
        if (err < tol * scale) { st = QD_OK; break; }
    }

    if (st != QD_OK) {
        double rj = qd_maxreljump_m(fv, N);
        st = (rj > QD_REL_JUMP_BRANCHCUT || decay_min > QD_DECAY_BRANCHCUT)
                 ? QD_BRANCHCUT : QD_NOCONV;
    }

    free(fv);
    R.value = value; R.abs_err = abs_err; R.status = st;
    R.n_used = N; R.radius = r;
    return R;
}

/* Adaptive-radius search: Fornberg-style downhill walk on the penalty,
 * starting at r = 1 and bracketing with a shrinking step ratio. */
static QdResultMachine qd_auto_m(QdSampleMachine f, void* ctx,
                                 double _Complex z0,
                                 double prec_goal_digits, int max_recursion) {
    double r = 1.0, step = 1.6;
    QdResultMachine cur = qd_fixed_m(f, ctx, z0, r, prec_goal_digits, max_recursion);
    double cur_pen = qd_penalty(cur.status, cur.abs_err, cur.n_used);
    QdResultMachine best = cur;
    double best_pen = cur_pen;

    for (int iter = 0; iter < 12 && step > 1.05; iter++) {
        QdResultMachine up = qd_fixed_m(f, ctx, z0, r * step, prec_goal_digits, max_recursion);
        QdResultMachine dn = qd_fixed_m(f, ctx, z0, r / step, prec_goal_digits, max_recursion);
        double up_pen = qd_penalty(up.status, up.abs_err, up.n_used);
        double dn_pen = qd_penalty(dn.status, dn.abs_err, dn.n_used);

        if (up_pen < cur_pen && up_pen <= dn_pen) { r *= step; cur = up; cur_pen = up_pen; }
        else if (dn_pen < cur_pen)                { r /= step; cur = dn; cur_pen = dn_pen; }
        else { step = sqrt(step); continue; }

        if (cur_pen < best_pen) { best = cur; best_pen = cur_pen; }
    }
    return best;
}

QdResultMachine qd_contour_residue_machine(QdSampleMachine f, void* ctx,
                                           double _Complex z0,
                                           double radius,
                                           double prec_goal_digits,
                                           int max_recursion) {
    if (max_recursion < 0) max_recursion = 0;
    if (radius > 0.0)
        return qd_fixed_m(f, ctx, z0, radius, prec_goal_digits, max_recursion);
    return qd_auto_m(f, ctx, z0, prec_goal_digits, max_recursion);
}

/* ================================================================== *
 *  MPFR path                                                          *
 * ================================================================== */
#ifdef USE_MPFR

#define QRND MPFR_RNDN

/* File-local complex toolkit: a pair of mpfr_t. Alias-safe — inputs are
 * read into temporaries before any output component is written. */
typedef struct { mpfr_t re, im; } qcx;

static void qcx_init(qcx* z, mpfr_prec_t p)  { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void qcx_clear(qcx* z)                { mpfr_clear(z->re); mpfr_clear(z->im); }
static void qcx_set(qcx* d, const qcx* s)    { mpfr_set(d->re, s->re, QRND); mpfr_set(d->im, s->im, QRND); }

static void qcx_sub(qcx* o, const qcx* a, const qcx* b) {
    mpfr_sub(o->re, a->re, b->re, QRND);
    mpfr_sub(o->im, a->im, b->im, QRND);
}

static void qcx_mul(qcx* o, const qcx* a, const qcx* b, mpfr_prec_t p) {
    mpfr_t ac, bd, ad, bc;
    mpfr_inits2(p, ac, bd, ad, bc, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, QRND);
    mpfr_mul(bd, a->im, b->im, QRND);
    mpfr_mul(ad, a->re, b->im, QRND);
    mpfr_mul(bc, a->im, b->re, QRND);
    mpfr_sub(o->re, ac, bd, QRND);
    mpfr_add(o->im, ad, bc, QRND);
    mpfr_clears(ac, bd, ad, bc, (mpfr_ptr)0);
}

/* o = a / b, Smith's algorithm. */
static void qcx_div(qcx* o, const qcx* a, const qcx* b, mpfr_prec_t p) {
    mpfr_t r, den, t1, t2;
    mpfr_inits2(p, r, den, t1, t2, (mpfr_ptr)0);
    if (mpfr_cmpabs(b->re, b->im) >= 0) {
        mpfr_div(r, b->im, b->re, QRND);          /* r = bi/br */
        mpfr_mul(t1, b->im, r, QRND);
        mpfr_add(den, b->re, t1, QRND);           /* den = br + bi r */
        mpfr_mul(t1, a->im, r, QRND);
        mpfr_add(t2, a->re, t1, QRND);
        mpfr_div(o->re, t2, den, QRND);           /* (ar + ai r)/den */
        mpfr_mul(t1, a->re, r, QRND);
        mpfr_sub(t2, a->im, t1, QRND);
        mpfr_div(o->im, t2, den, QRND);           /* (ai - ar r)/den */
    } else {
        mpfr_div(r, b->re, b->im, QRND);          /* r = br/bi */
        mpfr_mul(t1, b->re, r, QRND);
        mpfr_add(den, b->im, t1, QRND);           /* den = bi + br r */
        mpfr_mul(t1, a->re, r, QRND);
        mpfr_add(t2, t1, a->im, QRND);
        mpfr_div(o->re, t2, den, QRND);           /* (ar r + ai)/den */
        mpfr_mul(t1, a->im, r, QRND);
        mpfr_sub(t2, t1, a->re, QRND);
        mpfr_div(o->im, t2, den, QRND);           /* (ai r - ar)/den */
    }
    mpfr_clears(r, den, t1, t2, (mpfr_ptr)0);
}

static void qcx_abs(mpfr_t mag, const qcx* z) { mpfr_hypot(mag, z->re, z->im, QRND); }

/* Trapezoidal residue estimate into (out_re,out_im) (already init2'd at p).
 * fv[k] are the N samples; r the radius. */
static void qd_trap_mpfr(qcx* out, const qcx* fv, int N, double r, mpfr_prec_t p) {
    mpfr_t th, c, s, acc_re, acc_im, t;
    mpfr_inits2(p, th, c, s, acc_re, acc_im, t, (mpfr_ptr)0);
    mpfr_set_zero(acc_re, +1);
    mpfr_set_zero(acc_im, +1);
    mpfr_t two_pi; mpfr_init2(two_pi, p);
    mpfr_const_pi(two_pi, QRND);
    mpfr_mul_2ui(two_pi, two_pi, 1, QRND);        /* 2 pi */
    for (int k = 0; k < N; k++) {
        mpfr_mul_ui(th, two_pi, (unsigned long)k, QRND);
        mpfr_div_ui(th, th, (unsigned long)N, QRND);
        mpfr_sin_cos(s, c, th, QRND);
        /* acc += f_k * (c + i s) */
        mpfr_mul(t, fv[k].re, c, QRND);
        mpfr_fms(t, fv[k].im, s, t, QRND);        /* fi*s - fr*c */
        mpfr_sub(acc_re, acc_re, t, QRND);        /* += fr*c - fi*s */
        mpfr_mul(t, fv[k].re, s, QRND);
        mpfr_fma(t, fv[k].im, c, t, QRND);        /* fr*s + fi*c */
        mpfr_add(acc_im, acc_im, t, QRND);
    }
    /* out = acc * (r/N) */
    double rn = r / (double)N;
    mpfr_mul_d(out->re, acc_re, rn, QRND);
    mpfr_mul_d(out->im, acc_im, rn, QRND);
    mpfr_clears(th, c, s, acc_re, acc_im, t, two_pi, (mpfr_ptr)0);
}

/* Aitken Delta^2 on three complex estimates (MPFR). Writes the accelerated
 * value into `out`; sets *err to |out - s2| as a double gauge. */
static void qd_aitken_mpfr(qcx* out, const qcx* s0, const qcx* s1,
                           const qcx* s2, double* err, mpfr_prec_t p) {
    qcx d1, d2, den, q;
    qcx_init(&d1, p); qcx_init(&d2, p); qcx_init(&den, p); qcx_init(&q, p);
    qcx_sub(&d1, s1, s0);
    qcx_sub(&d2, s2, s1);
    qcx_sub(&den, &d2, &d1);
    mpfr_t aden, ad1, ad2, thresh;
    mpfr_inits2(p, aden, ad1, ad2, thresh, (mpfr_ptr)0);
    qcx_abs(aden, &den);
    qcx_abs(ad1, &d1);
    qcx_abs(ad2, &d2);
    mpfr_add(thresh, ad1, ad2, QRND);
    mpfr_mul_d(thresh, thresh, 1e-12, QRND);
    if (mpfr_zero_p(aden) || mpfr_cmp(aden, thresh) < 0) {
        qcx_set(out, s2);
        *err = mpfr_get_d(ad2, QRND);
    } else {
        qcx_mul(&q, &d2, &d2, p);     /* d2^2 */
        qcx_div(&q, &q, &den, p);     /* d2^2 / den */
        qcx_sub(out, s2, &q);         /* s2 - d2^2/den */
        mpfr_t e; mpfr_init2(e, p);
        qcx_sub(&q, out, s2);
        qcx_abs(e, &q);
        *err = mpfr_get_d(e, QRND);
        mpfr_clear(e);
    }
    mpfr_clears(aden, ad1, ad2, thresh, (mpfr_ptr)0);
    qcx_clear(&d1); qcx_clear(&d2); qcx_clear(&den); qcx_clear(&q);
}

/* Sample f at node `idx` of `N` on the circle of radius r about z0. */
static bool qd_sample_node_mpfr(QdSampleMpfr f, void* ctx,
                                const mpfr_t z0_re, const mpfr_t z0_im,
                                double r, int idx, int N,
                                qcx* out, mpfr_prec_t p) {
    mpfr_t th, c, s, zre, zim, two_pi;
    mpfr_inits2(p, th, c, s, zre, zim, two_pi, (mpfr_ptr)0);
    mpfr_const_pi(two_pi, QRND);
    mpfr_mul_2ui(two_pi, two_pi, 1, QRND);
    mpfr_mul_ui(th, two_pi, (unsigned long)idx, QRND);
    mpfr_div_ui(th, th, (unsigned long)N, QRND);
    mpfr_sin_cos(s, c, th, QRND);
    mpfr_mul_d(c, c, r, QRND);
    mpfr_mul_d(s, s, r, QRND);
    mpfr_add(zre, z0_re, c, QRND);
    mpfr_add(zim, z0_im, s, QRND);
    bool ok = f(ctx, zre, zim, out->re, out->im);
    if (ok && (!mpfr_number_p(out->re) || !mpfr_number_p(out->im))) ok = false;
    mpfr_clears(th, c, s, zre, zim, two_pi, (mpfr_ptr)0);
    return ok;
}

/* Relative-jump discontinuity gauge for the MPFR ring (computed via double
 * magnitudes — a coarse gauge is all that is needed). */
static double qd_maxreljump_mpfr(const qcx* f, int N) {
    double m = 0.0;
    for (int k = 0; k < N; k++) {
        const qcx* a = &f[k];
        const qcx* b = &f[(k + 1) % N];
        double ar = mpfr_get_d(a->re, QRND), ai = mpfr_get_d(a->im, QRND);
        double br = mpfr_get_d(b->re, QRND), bi = mpfr_get_d(b->im, QRND);
        double den = hypot(ar, ai) + hypot(br, bi);
        if (den > 0.0 && isfinite(den)) {
            double rj = hypot(br - ar, bi - ai) / den;
            if (rj > m) m = rj;
        }
    }
    return m;
}

static QdResultMpfr qd_fixed_mpfr(QdSampleMpfr f, void* ctx,
                                  const mpfr_t z0_re, const mpfr_t z0_im,
                                  double r, mpfr_prec_t p,
                                  double prec_goal_digits, int max_recursion,
                                  mpfr_t out_re, mpfr_t out_im) {
    QdResultMpfr R;
    R.status = QD_NONNUMERIC; R.abs_err = INFINITY; R.n_used = 0; R.radius = r;

    /* Initial N grows with precision: N0 ~ ln(10)*digits / 0.5, rounded up
     * to a power of two (>= the machine default). */
    int N0 = QD_N0_MACHINE;
    {
        double want = 4.61 * prec_goal_digits;   /* ln(10)/0.5 ~ 4.605 */
        int n = QD_N0_MACHINE;
        while ((double)n < want) n *= 2;
        N0 = n;
    }
    int N = N0;

    qcx* fv = malloc(sizeof(qcx) * (size_t)N);
    if (!fv) return R;
    for (int k = 0; k < N; k++) {
        qcx_init(&fv[k], p);
        if (!qd_sample_node_mpfr(f, ctx, z0_re, z0_im, r, k, N, &fv[k], p)) {
            for (int j = 0; j <= k; j++) qcx_clear(&fv[j]);
            free(fv);
            return R;
        }
    }

    qcx s_a, s_b, s_c, value, est;
    qcx_init(&s_a, p); qcx_init(&s_b, p); qcx_init(&s_c, p);
    qcx_init(&value, p); qcx_init(&est, p);
    qd_trap_mpfr(&s_c, fv, N, r, p);
    qcx_set(&value, &s_c);
    int have = 1;
    double prev_err = INFINITY, decay_min = 1.0, abs_err = INFINITY;
    QdStatus st = QD_NOCONV;
    double tol = pow(10.0, -prec_goal_digits);

    mpfr_t mag_diff, mag_c; mpfr_inits2(p, mag_diff, mag_c, (mpfr_ptr)0);

    for (int it = 0; it < max_recursion; it++) {
        int newN = 2 * N;
        /* Build the doubled sample set: evens are deep copies of the old
         * samples (cheap vs a function evaluation), odds are freshly
         * sampled. Then release the old array. */
        qcx* nv = malloc(sizeof(qcx) * (size_t)newN);
        if (!nv) break;
        bool ok = true;
        for (int i = 0; i < newN; i++) {
            qcx_init(&nv[i], p);
            if ((i & 1) == 0) {
                qcx_set(&nv[i], &fv[i / 2]);
            } else if (!qd_sample_node_mpfr(f, ctx, z0_re, z0_im, r, i, newN, &nv[i], p)) {
                ok = false;
                for (int j = 0; j <= i; j++) qcx_clear(&nv[j]);
                free(nv);
                break;
            }
        }
        if (!ok) { st = QD_NONNUMERIC; break; }
        for (int k = 0; k < N; k++) qcx_clear(&fv[k]);
        free(fv);
        fv = nv; N = newN;

        qcx S; qcx_init(&S, p);
        qd_trap_mpfr(&S, fv, N, r, p);
        qcx d; qcx_init(&d, p);
        qcx_sub(&d, &S, &s_c);
        qcx_abs(mag_diff, &d);
        double err = mpfr_get_d(mag_diff, QRND);
        qcx_clear(&d);

        double ratio = (isfinite(prev_err) && prev_err > 0.0) ? err / prev_err
                                                              : (err > 0.0 ? 1.0 : 0.0);
        if (ratio < decay_min) decay_min = ratio;
        prev_err = err;

        qcx_set(&s_a, &s_b); qcx_set(&s_b, &s_c); qcx_set(&s_c, &S);
        if (have < 3) have++;

        qcx_set(&est, &S);
        double est_err = err;
        if (have == 3) {
            qcx acc; qcx_init(&acc, p);
            double ae;
            qd_aitken_mpfr(&acc, &s_a, &s_b, &s_c, &ae, p);
            if (ae < est_err) { qcx_set(&est, &acc); est_err = ae; }
            qcx_clear(&acc);
        }
        qcx_set(&value, &est);
        abs_err = est_err;

        qcx_abs(mag_c, &s_c);
        double scale = 1.0 + mpfr_get_d(mag_c, QRND);
        qcx_clear(&S);
        /* Gate on the raw doubling error (see the machine-path note). */
        if (err < tol * scale) { st = QD_OK; break; }
    }

    if (st != QD_NONNUMERIC) {
        if (st != QD_OK) {
            double rj = qd_maxreljump_mpfr(fv, N);
            st = (rj > QD_REL_JUMP_BRANCHCUT || decay_min > QD_DECAY_BRANCHCUT)
                     ? QD_BRANCHCUT : QD_NOCONV;
        }
        mpfr_set(out_re, value.re, QRND);
        mpfr_set(out_im, value.im, QRND);
    }
    /* `fv` (the current sample array of N entries) is always live here:
     * the doubling-failure path frees only the half-built `nv`. */
    for (int k = 0; k < N; k++) qcx_clear(&fv[k]);
    free(fv);

    mpfr_clears(mag_diff, mag_c, (mpfr_ptr)0);
    qcx_clear(&s_a); qcx_clear(&s_b); qcx_clear(&s_c);
    qcx_clear(&value); qcx_clear(&est);

    R.status = st; R.abs_err = abs_err; R.n_used = N; R.radius = r;
    return R;
}

QdResultMpfr qd_contour_residue_mpfr(QdSampleMpfr f, void* ctx,
                                     const mpfr_t z0_re, const mpfr_t z0_im,
                                     double radius, long bits,
                                     double prec_goal_digits,
                                     int max_recursion,
                                     mpfr_t out_re, mpfr_t out_im) {
    if (max_recursion < 0) max_recursion = 0;
    mpfr_prec_t p = (mpfr_prec_t)bits;

    if (radius > 0.0)
        return qd_fixed_mpfr(f, ctx, z0_re, z0_im, radius, p,
                             prec_goal_digits, max_recursion, out_re, out_im);

    /* Adaptive radius: downhill walk on the penalty, as in the machine path,
     * but committing the best run's result into out_re/out_im. */
    mpfr_t br, bi; mpfr_init2(br, p); mpfr_init2(bi, p);
    double r = 1.0, step = 1.6;
    QdResultMpfr cur = qd_fixed_mpfr(f, ctx, z0_re, z0_im, r, p,
                                     prec_goal_digits, max_recursion, out_re, out_im);
    double cur_pen = qd_penalty(cur.status, cur.abs_err, cur.n_used);
    QdResultMpfr best = cur; double best_pen = cur_pen;
    mpfr_set(br, out_re, QRND); mpfr_set(bi, out_im, QRND);

    for (int iter = 0; iter < 12 && step > 1.05; iter++) {
        mpfr_t ur, ui, dr, di;
        mpfr_init2(ur, p); mpfr_init2(ui, p); mpfr_init2(dr, p); mpfr_init2(di, p);
        QdResultMpfr up = qd_fixed_mpfr(f, ctx, z0_re, z0_im, r * step, p,
                                        prec_goal_digits, max_recursion, ur, ui);
        QdResultMpfr dn = qd_fixed_mpfr(f, ctx, z0_re, z0_im, r / step, p,
                                        prec_goal_digits, max_recursion, dr, di);
        double up_pen = qd_penalty(up.status, up.abs_err, up.n_used);
        double dn_pen = qd_penalty(dn.status, dn.abs_err, dn.n_used);

        bool moved = true;
        if (up_pen < cur_pen && up_pen <= dn_pen) {
            r *= step; cur = up; cur_pen = up_pen;
            mpfr_set(out_re, ur, QRND); mpfr_set(out_im, ui, QRND);
        } else if (dn_pen < cur_pen) {
            r /= step; cur = dn; cur_pen = dn_pen;
            mpfr_set(out_re, dr, QRND); mpfr_set(out_im, di, QRND);
        } else { step = sqrt(step); moved = false; }

        if (moved && cur_pen < best_pen) {
            best = cur; best_pen = cur_pen;
            mpfr_set(br, out_re, QRND); mpfr_set(bi, out_im, QRND);
        }
        mpfr_clears(ur, ui, dr, di, (mpfr_ptr)0);
    }

    mpfr_set(out_re, br, QRND); mpfr_set(out_im, bi, QRND);
    mpfr_clears(br, bi, (mpfr_ptr)0);
    return best;
}

#endif /* USE_MPFR */

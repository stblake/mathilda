/*
 * matfunc_machine.c -- machine-precision kernels for the matrix functions.
 *
 *   mf_expm      exp(A): Pade scaling and squaring, degree 3/5/7/9/13 chosen by
 *                the Al-Mohy--Higham (2009) backward-error bounds, including the
 *                ell(A, m) correction that stops overscaling of nonnormal A.
 *   mf_expmv     exp(A) b: the Al-Mohy--Higham (2011) truncated Taylor action
 *                with the theta_m table and the trace shift; no exp(A) formed
 *                (it falls back to mf_expm . b when that is cheaper).
 *   mf_schur_fun the principal f(A) on the complex Schur form A = Q T Q^H:
 *                - log: inverse scaling and squaring (Al-Mohy--Higham 2012):
 *                  triangular square roots until T is near I, a Gauss--Legendre
 *                  partial-fraction Pade approximant of log(I + X), then the
 *                  diagonal and superdiagonal recomputed from closed forms;
 *                - sqrt: the Bjorck--Hammarling triangular recurrence;
 *                - any other f: Schur--Parlett (Davies--Higham 2003): the
 *                  eigenvalues clustered (delta = 0.1), the Schur form reordered
 *                  so each cluster is a contiguous block, each diagonal block by
 *                  a Taylor series about its mean eigenvalue, the off-diagonal
 *                  blocks by the block Parlett recurrence (triangular Sylvester
 *                  equations).
 *                A real matrix goes through the REAL Schur form (dgees) and is
 *                then triangularised (rsf2csf), so a real eigenvalue is exactly
 *                real: Log[-1.] lands on +I Pi, never on a rounding-decided
 *                branch.  A real answer is returned as real when every
 *                eigenvalue's f-value is consistent with a real matrix function.
 *
 * Every kernel BORROWS its inputs and writes caller-owned output.  The 1 x 1
 * case of every kernel is the scalar function itself (exp, log, z^s, f), so the
 * scalar REPL path, the packed / NDArray paths and Compile[] all reach the same
 * arithmetic.
 */

#include "matfunc_internal.h"
#include "lapack.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MF_U 1.1102230246251565e-16   /* unit roundoff 2^-53 */

/* ===================================================================== *
 *  Scalar complex arithmetic (no <complex.h>, per SPEC.md section 10).   *
 * ===================================================================== */

static inline void c_mul(double ar, double ai, double br, double bi,
                         double* r, double* i) {
    double t = ar * br - ai * bi;
    *i = ar * bi + ai * br;
    *r = t;
}

/* Smith's division. */
static inline void c_div(double ar, double ai, double br, double bi,
                         double* r, double* i) {
    if (fabs(br) >= fabs(bi)) {
        if (br == 0.0 && bi == 0.0) { *r = NAN; *i = NAN; return; }
        double q = bi / br, d = br + bi * q;
        double tr = (ar + ai * q) / d;
        *i = (ai - ar * q) / d;
        *r = tr;
    } else {
        double q = br / bi, d = bi + br * q;
        double tr = (ar * q + ai) / d;
        *i = (ai * q - ar) / d;
        *r = tr;
    }
}

static inline double c_abs(double r, double i) { return hypot(r, i); }

/* Principal log; a -0. imaginary part is read as +0., so a negative real
 * argument always lands on +I Pi (the reference's Log[-1.] == I Pi). */
static void c_log(double zr, double zi, double* r, double* i) {
    if (zi == 0.0) zi = 0.0;
    *r = log(hypot(zr, zi));
    *i = atan2(zi, zr);
}

static void c_exp(double zr, double zi, double* r, double* i) {
    double e = exp(zr);
    if (zi == 0.0) { *r = e; *i = 0.0; return; }
    *r = e * cos(zi);
    *i = e * sin(zi);
}

/* Principal square root. */
static void c_sqrt(double zr, double zi, double* r, double* i) {
    if (zi == 0.0) zi = 0.0;
    if (zr == 0.0 && zi == 0.0) { *r = 0.0; *i = 0.0; return; }
    double t = sqrt((fabs(zr) + hypot(zr, zi)) * 0.5);
    if (zr >= 0.0) { *r = t; *i = zi / (2.0 * t); }
    else           { *r = fabs(zi) / (2.0 * t); *i = copysign(t, zi); }
}

/* log1p for a complex argument: log(1 + w) without cancelling when w is small. */
static void c_log1p(double wr, double wi, double* r, double* i) {
    *r = 0.5 * log1p(2.0 * wr + wr * wr + wi * wi);
    *i = atan2(wi, 1.0 + wr);
}

/* atanh(z) = (log1p(z) - log1p(-z)) / 2. */
static void c_atanh(double zr, double zi, double* r, double* i) {
    double ar, ai, br, bi;
    c_log1p(zr, zi, &ar, &ai);
    c_log1p(-zr, -zi, &br, &bi);
    *r = 0.5 * (ar - br);
    *i = 0.5 * (ai - bi);
}

/* ===================================================================== *
 *  The closed-form scalar functions (and their derivatives).            *
 * ===================================================================== */

int mf_fn_exp(void* ctx, double zr, double zi, int k, double* fr, double* fi) {
    (void)ctx; (void)k;
    c_exp(zr, zi, fr, fi);
    return (isfinite(*fr) && isfinite(*fi)) ? MF_OK : MF_FAIL;
}

int mf_fn_log(void* ctx, double zr, double zi, int k, double* fr, double* fi) {
    (void)ctx;
    if (zr == 0.0 && zi == 0.0) return MF_ZEROEIG;
    if (k == 0) { c_log(zr, zi, fr, fi); return MF_OK; }
    /* (-1)^(k-1) (k-1)! / z^k */
    double pr = 1.0, pi = 0.0;
    for (int j = 0; j < k; j++) c_mul(pr, pi, zr, zi, &pr, &pi);
    double fact = 1.0;
    for (int j = 2; j < k; j++) fact *= (double)j;
    if (k % 2 == 0) fact = -fact;
    c_div(fact, 0.0, pr, pi, fr, fi);
    return (isfinite(*fr) && isfinite(*fi)) ? MF_OK : MF_FAIL;
}

/* z^s, principal branch; ctx = double[2] {Re s, Im s}.  A positive real base
 * with a real exponent uses pow() (exactly the scalar Power), a negative or
 * complex base exp(s log z) -- (-4.)^0.5 is 1.2246467991473532*^-16 + 2. I,
 * as the scalar Power gives it. */
static void c_pow(double zr, double zi, double sr, double si, double* r, double* i) {
    if (zi == 0.0 && si == 0.0 && zr > 0.0) { *r = pow(zr, sr); *i = 0.0; return; }
    double lr, li, er, ei;
    c_log(zr, zi, &lr, &li);
    c_mul(sr, si, lr, li, &er, &ei);
    c_exp(er, ei, r, i);
}

int mf_fn_pow(void* ctx, double zr, double zi, int k, double* fr, double* fi) {
    const double* s = (const double*)ctx;
    /* coefficient s (s - 1) ... (s - k + 1) */
    double cr = 1.0, ci = 0.0;
    for (int j = 0; j < k; j++) c_mul(cr, ci, s[0] - (double)j, s[1], &cr, &ci);
    if (cr == 0.0 && ci == 0.0) { *fr = 0.0; *fi = 0.0; return MF_OK; }
    if (zr == 0.0 && zi == 0.0) {
        double er = s[0] - (double)k;
        if (er > 0.0) { *fr = 0.0; *fi = 0.0; return MF_OK; }
        if (er == 0.0 && s[1] == 0.0) { *fr = cr; *fi = ci; return MF_OK; }
        return (k == 0) ? MF_ZEROEIG : MF_NOSOL;
    }
    double pr, pi;
    c_pow(zr, zi, s[0] - (double)k, s[1], &pr, &pi);
    c_mul(cr, ci, pr, pi, fr, fi);
    return (isfinite(*fr) && isfinite(*fi)) ? MF_OK : MF_FAIL;
}

/* ===================================================================== *
 *  Dense helpers over real (w = 1) or interleaved complex (w = 2),      *
 *  column-major n x n buffers.                                          *
 * ===================================================================== */

/* C = A B. */
static void mm(int n, bool cx, const double* A, const double* B, double* C) {
#ifdef USE_LAPACK
    if (cx) {
        double one[2] = { 1.0, 0.0 }, zero[2] = { 0.0, 0.0 };
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
                    (void*)one, (const void*)A, n, (const void*)B, n,
                    (void*)zero, (void*)C, n);
    } else {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
                    1.0, A, n, B, n, 0.0, C, n);
    }
#else
    size_t nn = (size_t)n * (size_t)n;
    memset(C, 0, (cx ? 2 : 1) * nn * sizeof(double));
    for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++) {
            size_t bk = (size_t)k + (size_t)j * (size_t)n;
            if (cx) {
                double br = B[2 * bk], bi = B[2 * bk + 1];
                if (br == 0.0 && bi == 0.0) continue;
                for (int i = 0; i < n; i++) {
                    size_t ai = (size_t)i + (size_t)k * (size_t)n;
                    size_t ci = (size_t)i + (size_t)j * (size_t)n;
                    double ar = A[2 * ai], aim = A[2 * ai + 1];
                    C[2 * ci]     += ar * br - aim * bi;
                    C[2 * ci + 1] += ar * bi + aim * br;
                }
            } else {
                double b = B[bk];
                if (b == 0.0) continue;
                for (int i = 0; i < n; i++)
                    C[(size_t)i + (size_t)j * (size_t)n] += A[(size_t)i + (size_t)k * (size_t)n] * b;
            }
        }
#endif
}

/* y = A x (complex when cx; x, y the same element type). */
static void mv(int n, bool cx, const double* A, const double* x, double* y) {
#ifdef USE_LAPACK
    if (cx) {
        double one[2] = { 1.0, 0.0 }, zero[2] = { 0.0, 0.0 };
        cblas_zgemv(CblasColMajor, CblasNoTrans, n, n, (void*)one, (const void*)A, n,
                    (const void*)x, 1, (void*)zero, (void*)y, 1);
    } else {
        cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, A, n, x, 1, 0.0, y, 1);
    }
#else
    memset(y, 0, (cx ? 2 : 1) * (size_t)n * sizeof(double));
    for (int k = 0; k < n; k++)
        for (int i = 0; i < n; i++) {
            size_t a = (size_t)i + (size_t)k * (size_t)n;
            if (cx) {
                y[2 * i]     += A[2 * a] * x[2 * k] - A[2 * a + 1] * x[2 * k + 1];
                y[2 * i + 1] += A[2 * a] * x[2 * k + 1] + A[2 * a + 1] * x[2 * k];
            } else {
                y[i] += A[a] * x[k];
            }
        }
#endif
}

/* Induced 1-norm: the largest column sum of moduli. */
static double norm1(int n, bool cx, const double* A) {
    double best = 0.0;
    for (int j = 0; j < n; j++) {
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            size_t k = (size_t)i + (size_t)j * (size_t)n;
            s += cx ? c_abs(A[2 * k], A[2 * k + 1]) : fabs(A[k]);
        }
        if (s > best || isnan(s)) best = s;
    }
    return best;
}

static double norm_inf_vec(int n, bool cx, const double* x) {
    double best = 0.0;
    for (int i = 0; i < n; i++) {
        double a = cx ? c_abs(x[2 * i], x[2 * i + 1]) : fabs(x[i]);
        if (a > best) best = a;
    }
    return best;
}

static bool all_finite(size_t count, const double* A) {
    for (size_t k = 0; k < count; k++) if (!isfinite(A[k])) return false;
    return true;
}

/* Solve P X = Q in place (Q becomes X); P is destroyed.  LAPACK when present,
 * else Gaussian elimination with partial pivoting.  0 on success. */
static int solve(int n, bool cx, double* P, double* Q) {
#ifdef USE_LAPACK
    if (mathilda_lapack_probe()) {
        int* ipiv = (int*)malloc(sizeof(int) * (size_t)n);
        if (!ipiv) return -1;
        int info = cx ? mat_lapack_zgesv(n, n, P, n, ipiv, Q, n)
                      : mat_lapack_dgesv(n, n, P, n, ipiv, Q, n);
        free(ipiv);
        return info;
    }
#endif
    int w = cx ? 2 : 1;
#define PE(i, j) (P + (size_t)w * ((size_t)(i) + (size_t)(j) * (size_t)n))
#define QE(i, j) (Q + (size_t)w * ((size_t)(i) + (size_t)(j) * (size_t)n))
    for (int k = 0; k < n; k++) {
        int piv = k;
        double best = -1.0;
        for (int i = k; i < n; i++) {
            double a = cx ? c_abs(PE(i, k)[0], PE(i, k)[1]) : fabs(PE(i, k)[0]);
            if (a > best) { best = a; piv = i; }
        }
        if (best == 0.0) return 1;
        if (piv != k) {
            for (int j = 0; j < n; j++) {
                for (int c = 0; c < w; c++) {
                    double t = PE(k, j)[c]; PE(k, j)[c] = PE(piv, j)[c]; PE(piv, j)[c] = t;
                    t = QE(k, j)[c]; QE(k, j)[c] = QE(piv, j)[c]; QE(piv, j)[c] = t;
                }
            }
        }
        for (int i = k + 1; i < n; i++) {
            double lr, li = 0.0;
            if (cx) c_div(PE(i, k)[0], PE(i, k)[1], PE(k, k)[0], PE(k, k)[1], &lr, &li);
            else    lr = PE(i, k)[0] / PE(k, k)[0];
            if (lr == 0.0 && li == 0.0) continue;
            for (int j = k; j < n; j++) {
                if (cx) {
                    double pr, pi;
                    c_mul(lr, li, PE(k, j)[0], PE(k, j)[1], &pr, &pi);
                    PE(i, j)[0] -= pr; PE(i, j)[1] -= pi;
                } else PE(i, j)[0] -= lr * PE(k, j)[0];
            }
            for (int j = 0; j < n; j++) {
                if (cx) {
                    double pr, pi;
                    c_mul(lr, li, QE(k, j)[0], QE(k, j)[1], &pr, &pi);
                    QE(i, j)[0] -= pr; QE(i, j)[1] -= pi;
                } else QE(i, j)[0] -= lr * QE(k, j)[0];
            }
        }
    }
    for (int j = 0; j < n; j++)
        for (int i = n - 1; i >= 0; i--) {
            double sr = QE(i, j)[0], si = cx ? QE(i, j)[1] : 0.0;
            for (int k = i + 1; k < n; k++) {
                if (cx) {
                    double pr, pi;
                    c_mul(PE(i, k)[0], PE(i, k)[1], QE(k, j)[0], QE(k, j)[1], &pr, &pi);
                    sr -= pr; si -= pi;
                } else sr -= PE(i, k)[0] * QE(k, j)[0];
            }
            if (cx) c_div(sr, si, PE(i, i)[0], PE(i, i)[1], &QE(i, j)[0], &QE(i, j)[1]);
            else    QE(i, j)[0] = sr / PE(i, i)[0];
        }
#undef PE
#undef QE
    return 0;
}

/* ===================================================================== *
 *  exp(A): Pade scaling and squaring (Al-Mohy & Higham 2009).           *
 * ===================================================================== */

static const double PADE3[4]  = { 120., 60., 12., 1. };
static const double PADE5[6]  = { 30240., 15120., 3360., 420., 30., 1. };
static const double PADE7[8]  = { 17297280., 8648640., 1995840., 277200., 25200.,
                                  1512., 56., 1. };
static const double PADE9[10] = { 17643225600., 8821612800., 2075673600., 302702400.,
                                  30270240., 2162160., 110880., 3960., 90., 1. };
static const double PADE13[14] = { 64764752532480000., 32382376266240000.,
                                   7771770303897600., 1187353796428800.,
                                   129060195264000., 10559470521600., 670442572800.,
                                   33522128640., 1323241920., 40840800., 960960.,
                                   16380., 182., 1. };

/* 1-norm of |A|^p for the nonnegative |A| (exact: 1^T |A|^p by p row-vector
 * products), times scale^p. */
static double abs_power_norm(int n, const double* absA, int p) {
    double* v = (double*)malloc(sizeof(double) * (size_t)n);
    double* t = (double*)malloc(sizeof(double) * (size_t)n);
    if (!v || !t) { free(v); free(t); return INFINITY; }
    for (int i = 0; i < n; i++) v[i] = 1.0;
    for (int q = 0; q < p; q++) {
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int i = 0; i < n; i++) s += v[i] * absA[(size_t)i + (size_t)j * (size_t)n];
            t[j] = s;
        }
        memcpy(v, t, sizeof(double) * (size_t)n);
    }
    double best = 0.0;
    for (int i = 0; i < n; i++) if (v[i] > best) best = v[i];
    free(v); free(t);
    return best;
}

/* The ell(A, m) correction of Al-Mohy & Higham (2009, eq. 3.4 ff.): the
 * number of extra squarings the degree-m approximant needs for A scaled by
 * 2^-s.  `nA` is ||A||_1 of the unscaled A. */
static int ell(int n, const double* absA, double nA, int m, int s) {
    int p = 2 * m + 1;
    /* |c_{2m+1}|^-1 = binomial(2m, m) * (2m+1)! */
    double binom = 1.0;
    for (int j = 1; j <= m; j++) binom = binom * (double)(m + j) / (double)j;
    double fact = 1.0;
    for (int j = 2; j <= p; j++) fact *= (double)j;
    double abs_c_recip = binom * fact;
    double nrm = abs_power_norm(n, absA, p);
    if (!(nrm > 0.0)) return 0;
    /* alpha for 2^-s A: 2^-sp ||A^p|| / (2^-s ||A|| c^-1) */
    double alpha = ldexp(nrm / (nA * abs_c_recip), -s * (p - 1));
    if (!(alpha > 0.0)) return 0;
    double v = ceil(log2(alpha / MF_U) / (2.0 * m));
    return v > 0.0 ? (int)v : 0;
}

/* E = sum_j c[j] P[j] (+ c0 I), element-wise. */
static void combine(int n, bool cx, int count, const double* const* P,
                    const double* c, double c_identity, double* E) {
    size_t nn = (size_t)n * (size_t)n, sz = (cx ? 2 : 1) * nn;
    memset(E, 0, sz * sizeof(double));
    for (int j = 0; j < count; j++) {
        double cj = c[j];
        const double* Pj = P[j];
        for (size_t k = 0; k < sz; k++) E[k] += cj * Pj[k];
    }
    for (int i = 0; i < n; i++)
        E[(cx ? 2 : 1) * ((size_t)i + (size_t)i * (size_t)n)] += c_identity;
}

int mf_expm(const double* A, int n, bool cx, double* E) {
    int w = cx ? 2 : 1;
    size_t nn = (size_t)n * (size_t)n, sz = (size_t)w * nn;
    if (n < 1 || !all_finite(sz, A)) return MF_FAIL;

    if (n == 1) {                       /* the scalar exp itself */
        if (cx) c_exp(A[0], A[1], &E[0], &E[1]);
        else    E[0] = exp(A[0]);
        return all_finite(sz, E) ? MF_OK : MF_FAIL;
    }

    double nA = norm1(n, cx, A);
    if (nA == 0.0) {                    /* exp(0) = I */
        memset(E, 0, sz * sizeof(double));
        for (int i = 0; i < n; i++) E[(size_t)w * ((size_t)i + (size_t)i * (size_t)n)] = 1.0;
        return MF_OK;
    }

    double* buf = (double*)malloc(sizeof(double) * sz * 8 + sizeof(double) * nn);
    if (!buf) return MF_FAIL;
    double *A2 = buf, *A4 = buf + sz, *A6 = buf + 2 * sz, *A8 = buf + 3 * sz;
    double *U = buf + 4 * sz, *V = buf + 5 * sz, *T1 = buf + 6 * sz, *T2 = buf + 7 * sz;
    double* absA = buf + 8 * sz;
    for (size_t k = 0; k < nn; k++)
        absA[k] = cx ? c_abs(A[2 * k], A[2 * k + 1]) : fabs(A[k]);

    mm(n, cx, A, A, A2);
    mm(n, cx, A2, A2, A4);
    mm(n, cx, A4, A2, A6);
    double d4 = pow(norm1(n, cx, A4), 0.25);
    double d6 = pow(norm1(n, cx, A6), 1.0 / 6.0);
    double eta1 = fmax(d4, d6);

    const double* pw[5] = { NULL, A2, A4, A6, A8 };
    const double* coef = NULL;
    int deg = 0;
    if (eta1 <= 1.495585217958292e-2 && ell(n, absA, nA, 3, 0) == 0) { deg = 3; coef = PADE3; }
    else if (eta1 <= 2.539398330063230e-1 && ell(n, absA, nA, 5, 0) == 0) { deg = 5; coef = PADE5; }
    double d8 = 0.0, eta3 = 0.0;
    if (!deg) {
        mm(n, cx, A4, A4, A8);
        d8 = pow(norm1(n, cx, A8), 0.125);
        eta3 = fmax(d6, d8);
        if (eta3 <= 9.504178996162932e-1 && ell(n, absA, nA, 7, 0) == 0) { deg = 7; coef = PADE7; }
        else if (eta3 <= 2.097847961257068e0 && ell(n, absA, nA, 9, 0) == 0) { deg = 9; coef = PADE9; }
    }

    int s = 0;
    int status = MF_OK;
    if (deg) {
        /* U = A (sum_j b_{2j+1} A^{2j}), V = sum_j b_{2j} A^{2j}. */
        int half = (deg - 1) / 2;       /* number of even powers past I */
        double codd[5], ceven[5];
        for (int j = 1; j <= half; j++) { codd[j - 1] = coef[2 * j + 1]; ceven[j - 1] = coef[2 * j]; }
        combine(n, cx, half, pw + 1, codd, coef[1], T1);
        mm(n, cx, A, T1, U);
        combine(n, cx, half, pw + 1, ceven, coef[0], V);
    } else {
        mm(n, cx, A4, A6, T1);
        double d10 = pow(norm1(n, cx, T1), 0.1);
        double eta4 = fmax(d8, d10);
        double eta5 = fmin(eta3, eta4);
        if (eta5 > 0.0) {
            double sd = ceil(log2(eta5 / 4.25));
            if (!isfinite(sd) || sd > 1100.0) { free(buf); return MF_FAIL; }
            s = sd > 0.0 ? (int)sd : 0;
        }
        s += ell(n, absA, nA, 13, s);
        /* Scale A, A^2, A^4, A^6 by 2^-s, 2^-2s, 2^-4s, 2^-6s. */
        double* As = A8;                /* A8 is free in the degree-13 branch */
        for (size_t k = 0; k < sz; k++) {
            As[k] = ldexp(A[k], -s);
            A2[k] = ldexp(A2[k], -2 * s);
            A4[k] = ldexp(A4[k], -4 * s);
            A6[k] = ldexp(A6[k], -6 * s);
        }
        const double* P3[3] = { A6, A4, A2 };
        /* U = As (A6 (b13 A6 + b11 A4 + b9 A2) + b7 A6 + b5 A4 + b3 A2 + b1 I) */
        double c1[3] = { PADE13[13], PADE13[11], PADE13[9] };
        combine(n, cx, 3, P3, c1, 0.0, T1);
        mm(n, cx, A6, T1, T2);
        double c2[3] = { PADE13[7], PADE13[5], PADE13[3] };
        combine(n, cx, 3, P3, c2, PADE13[1], T1);
        for (size_t k = 0; k < sz; k++) T1[k] += T2[k];
        mm(n, cx, As, T1, U);
        /* V = A6 (b12 A6 + b10 A4 + b8 A2) + b6 A6 + b4 A4 + b2 A2 + b0 I */
        double c3[3] = { PADE13[12], PADE13[10], PADE13[8] };
        combine(n, cx, 3, P3, c3, 0.0, T1);
        mm(n, cx, A6, T1, T2);
        double c4[3] = { PADE13[6], PADE13[4], PADE13[2] };
        combine(n, cx, 3, P3, c4, PADE13[0], V);
        for (size_t k = 0; k < sz; k++) V[k] += T2[k];
    }

    /* (V - U) X = (V + U) */
    for (size_t k = 0; k < sz; k++) {
        T1[k] = V[k] - U[k];
        E[k]  = V[k] + U[k];
    }
    if (solve(n, cx, T1, E) != 0) status = MF_FAIL;
    for (int q = 0; q < s && status == MF_OK; q++) {
        mm(n, cx, E, E, T2);
        memcpy(E, T2, sz * sizeof(double));
    }
    free(buf);
    if (status == MF_OK && !all_finite(sz, E)) status = MF_FAIL;
    return status;
}

/* ===================================================================== *
 *  exp(A) b: the truncated-Taylor action (Al-Mohy & Higham 2011).       *
 * ===================================================================== */

/* theta_m for the unit roundoff 2^-53: m = 1..30 from Higham's "Functions of
 * Matrices" table A.3, then 35..55 from Al-Mohy & Higham (2011) table 3.1. */
const int    EXPMV_M[]     = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                      16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
                                      29, 30, 35, 40, 45, 50, 55 };
const double EXPMV_THETA[] = { 2.29e-16, 2.58e-8, 1.39e-5, 3.40e-4, 2.40e-3, 9.07e-3,
                                      2.38e-2, 5.00e-2, 8.96e-2, 1.44e-1, 2.14e-1, 3.00e-1,
                                      4.00e-1, 5.14e-1, 6.41e-1, 7.81e-1, 9.31e-1, 1.09,
                                      1.26, 1.44, 1.62, 1.82, 2.01, 2.22, 2.43, 2.64, 2.86,
                                      3.08, 3.31, 3.54, 4.7, 6.0, 7.2, 8.5, 9.9 };

const size_t EXPMV_NTHETA = sizeof EXPMV_M / sizeof EXPMV_M[0];

int mf_expmv(const double* A, int n, bool cx, const double* b, bool bcx, double* out) {
    bool z = cx || bcx;
    int w = z ? 2 : 1;
    size_t nn = (size_t)n * (size_t)n;
    if (n < 1 || !all_finite((cx ? 2 : 1) * nn, A) || !all_finite((size_t)(bcx ? 2 : 1) * (size_t)n, b))
        return MF_FAIL;

    /* Promote both operands to the working element type. */
    double* M = (double*)malloc(sizeof(double) * (size_t)w * nn);
    double* x = (double*)malloc(sizeof(double) * (size_t)w * (size_t)n);
    double* y = (double*)malloc(sizeof(double) * (size_t)w * (size_t)n);
    if (!M || !x || !y) { free(M); free(x); free(y); return MF_FAIL; }
    for (size_t k = 0; k < nn; k++) {
        if (z) { M[2 * k] = cx ? A[2 * k] : A[k]; M[2 * k + 1] = cx ? A[2 * k + 1] : 0.0; }
        else M[k] = A[k];
    }
    for (int i = 0; i < n; i++) {
        if (z) { x[2 * i] = bcx ? b[2 * i] : b[i]; x[2 * i + 1] = bcx ? b[2 * i + 1] : 0.0; }
        else x[i] = b[i];
    }

    /* Shift by mu = trace(A)/n. */
    double mur = 0.0, mui = 0.0;
    for (int i = 0; i < n; i++) {
        size_t d = (size_t)i + (size_t)i * (size_t)n;
        mur += z ? M[2 * d] : M[d];
        if (z) mui += M[2 * d + 1];
    }
    mur /= n; mui /= n;
    for (int i = 0; i < n; i++) {
        size_t d = (size_t)i + (size_t)i * (size_t)n;
        if (z) { M[2 * d] -= mur; M[2 * d + 1] -= mui; } else M[d] -= mur;
    }
    double anorm = norm1(n, z, M);

    /* Choose (m, s) minimising the matvec count m s. */
    int mstar = 0;
    double sbest = 1.0, cost = 0.0;
    if (anorm > 0.0) {
        double best = INFINITY;
        for (size_t t = 0; t < sizeof EXPMV_M / sizeof EXPMV_M[0]; t++) {
            double sc = ceil(anorm / EXPMV_THETA[t]);
            double c = (double)EXPMV_M[t] * sc;
            if (c < best) { best = c; mstar = EXPMV_M[t]; sbest = sc; }
        }
        cost = best;
    }

    int status = MF_OK;
    /* A dense exp(A) costs ~ 20 n^3 flops against 2 n^2 per matvec: past that
     * break-even the full exponential is the cheaper (and equally accurate)
     * route. */
    if (cost > 10.0 * (double)n || !isfinite(cost)) {
        double* E = (double*)malloc(sizeof(double) * (size_t)w * nn);
        double* Az = M;   /* reload the unshifted operand */
        for (size_t k = 0; k < nn; k++) {
            if (z) { Az[2 * k] = cx ? A[2 * k] : A[k]; Az[2 * k + 1] = cx ? A[2 * k + 1] : 0.0; }
            else Az[k] = A[k];
        }
        if (!E || mf_expm(Az, n, z, E) != MF_OK) status = MF_FAIL;
        else mv(n, z, E, x, out);
        free(E);
        free(M); free(x); free(y);
        if (status == MF_OK && !all_finite((size_t)w * (size_t)n, out)) status = MF_FAIL;
        return status;
    }

    int s = (int)sbest;
    double etar, etai;
    c_exp(mur / s, mui / s, &etar, &etai);
    double* F = out;
    memcpy(F, x, sizeof(double) * (size_t)w * (size_t)n);
    for (int i = 0; i < s; i++) {
        double c1 = norm_inf_vec(n, z, x);
        for (int j = 1; j <= mstar; j++) {
            mv(n, z, M, x, y);
            double scale = 1.0 / ((double)s * (double)j);
            for (int k = 0; k < w * n; k++) x[k] = y[k] * scale;
            double c2 = norm_inf_vec(n, z, x);
            for (int k = 0; k < w * n; k++) F[k] += x[k];
            if (c1 + c2 <= MF_U * norm_inf_vec(n, z, F)) break;
            c1 = c2;
        }
        for (int k = 0; k < n; k++) {
            if (z) c_mul(etar, etai, F[2 * k], F[2 * k + 1], &F[2 * k], &F[2 * k + 1]);
            else F[k] *= etar;
        }
        memcpy(x, F, sizeof(double) * (size_t)w * (size_t)n);
    }
    free(M); free(x); free(y);
    if (!all_finite((size_t)w * (size_t)n, out)) status = MF_FAIL;
    return status;
}

/* ===================================================================== *
 *  Schur-form machinery (complex upper-triangular T, unitary Q).        *
 * ===================================================================== */

#define CE(M, i, j) ((M) + 2 * ((size_t)(i) + (size_t)(j) * (size_t)n))

/* Triangularise a real quasi-triangular Schur form held as complex (T, Z):
 * each 2 x 2 bump is removed by a Givens rotation (the rsf2csf method), so the
 * 1 x 1 real eigenvalues stay exactly real. */
static void rsf2csf(int n, double* T, double* Z) {
    for (int m = n - 1; m >= 1; m--) {
        double* c = CE(T, m, m - 1);
        double a_r = CE(T, m - 1, m - 1)[0], a_i = CE(T, m - 1, m - 1)[1];
        double d_r = CE(T, m, m)[0], d_i = CE(T, m, m)[1];
        if (c_abs(c[0], c[1]) > MF_U * 2.0 * (c_abs(a_r, a_i) + c_abs(d_r, d_i))) {
            double b_r = CE(T, m - 1, m)[0], b_i = CE(T, m - 1, m)[1];
            /* eigenvalue of the 2 x 2 block, shifted by -d */
            double tr_r = a_r + d_r, tr_i = a_i + d_i;
            double det_r, det_i, bc_r, bc_i;
            c_mul(a_r, a_i, d_r, d_i, &det_r, &det_i);
            c_mul(b_r, b_i, c[0], c[1], &bc_r, &bc_i);
            det_r -= bc_r; det_i -= bc_i;
            double h_r = 0.5 * tr_r, h_i = 0.5 * tr_i, q_r, q_i;
            c_mul(h_r, h_i, h_r, h_i, &q_r, &q_i);
            q_r -= det_r; q_i -= det_i;
            double s_r, s_i;
            c_sqrt(q_r, q_i, &s_r, &s_i);
            double mu_r = h_r + s_r - d_r, mu_i = h_i + s_i - d_i;
            double r = hypot(c_abs(mu_r, mu_i), c_abs(c[0], c[1]));
            double cc_r = mu_r / r, cc_i = mu_i / r;
            double sn = c[0] / r;            /* the subdiagonal entry is real */
            /* rows m-1, m, columns m-1..n-1:  G = [[conj(cc), sn], [-sn, cc]] */
            for (int j = m - 1; j < n; j++) {
                double* x = CE(T, m - 1, j);
                double* y = CE(T, m, j);
                double xr = x[0], xi = x[1], yr = y[0], yi = y[1], t_r, t_i;
                c_mul(cc_r, -cc_i, xr, xi, &t_r, &t_i);
                x[0] = t_r + sn * yr; x[1] = t_i + sn * yi;
                c_mul(cc_r, cc_i, yr, yi, &t_r, &t_i);
                y[0] = t_r - sn * xr; y[1] = t_i - sn * xi;
            }
            /* columns m-1, m by G^H = [[cc, -sn], [sn, conj(cc)]] */
            for (int pass = 0; pass < 2; pass++) {
                double* M = pass == 0 ? T : Z;
                int rows = pass == 0 ? m + 1 : n;
                for (int i = 0; i < rows; i++) {
                    double* p = CE(M, i, m - 1);
                    double* q = CE(M, i, m);
                    double pr = p[0], pi = p[1], qr = q[0], qi = q[1], t_r, t_i;
                    c_mul(pr, pi, cc_r, cc_i, &t_r, &t_i);
                    p[0] = t_r + qr * sn; p[1] = t_i + qi * sn;
                    c_mul(qr, qi, cc_r, -cc_i, &t_r, &t_i);
                    q[0] = t_r - pr * sn; q[1] = t_i - pi * sn;
                }
            }
        }
        c[0] = 0.0; c[1] = 0.0;
    }
}

/* Complex Schur form of A (complex col-major output T, Q).  A real input goes
 * through the real Schur form (exactly-real real eigenvalues). */
static int schur_complex(const double* A, int n, bool cx, double* T, double* Q) {
    if (!mathilda_lapack_probe()) return MF_FAIL;
    size_t nn = (size_t)n * (size_t)n;
    if (cx) {
        double* wv = (double*)malloc(sizeof(double) * 2 * (size_t)n);
        if (!wv) return MF_FAIL;
        memcpy(T, A, sizeof(double) * 2 * nn);
        int info = mat_lapack_zgees(n, T, n, wv, Q, n);
        free(wv);
        if (info != 0) return MF_FAIL;
    } else {
        double* Tr = (double*)malloc(sizeof(double) * nn);
        double* Qr = (double*)malloc(sizeof(double) * nn);
        double* wr = (double*)malloc(sizeof(double) * (size_t)n);
        double* wi = (double*)malloc(sizeof(double) * (size_t)n);
        if (!Tr || !Qr || !wr || !wi) { free(Tr); free(Qr); free(wr); free(wi); return MF_FAIL; }
        memcpy(Tr, A, sizeof(double) * nn);
        int info = mat_lapack_dgees(n, Tr, n, wr, wi, Qr, n);
        if (info == 0)
            for (size_t k = 0; k < nn; k++) {
                T[2 * k] = Tr[k]; T[2 * k + 1] = 0.0;
                Q[2 * k] = Qr[k]; Q[2 * k + 1] = 0.0;
            }
        free(Tr); free(Qr); free(wr); free(wi);
        if (info != 0) return MF_FAIL;
        rsf2csf(n, T, Q);
    }
    /* the strictly lower triangle is zero by construction */
    for (int j = 0; j < n; j++)
        for (int i = j + 1; i < n; i++) { CE(T, i, j)[0] = 0.0; CE(T, i, j)[1] = 0.0; }
    return MF_OK;
}

/* Swap the adjacent diagonal entries k, k+1 of the triangular T, updating Q
 * (the ztrexc step: one Givens rotation from zlartg). */
static void schur_swap(int n, double* T, double* Q, int k) {
    double t11r = CE(T, k, k)[0], t11i = CE(T, k, k)[1];
    double t22r = CE(T, k + 1, k + 1)[0], t22i = CE(T, k + 1, k + 1)[1];
    double fr = CE(T, k, k + 1)[0], fi = CE(T, k, k + 1)[1];
    double gr = t22r - t11r, gi = t22i - t11i;
    double cs, snr, sni;
    double fa = c_abs(fr, fi), ga = c_abs(gr, gi);
    if (ga == 0.0) { cs = 1.0; snr = 0.0; sni = 0.0; }
    else if (fa == 0.0) { cs = 0.0; snr = gr / ga; sni = -gi / ga; }
    else {
        double d = hypot(fa, ga);
        cs = fa / d;
        /* sn = (f / |f|) conj(g) / d */
        c_mul(fr / fa, fi / fa, gr / d, -gi / d, &snr, &sni);
    }
    /* rows k, k+1, columns k+2..n-1: x' = cs x + sn y, y' = cs y - conj(sn) x */
    for (int j = k + 2; j < n; j++) {
        double* x = CE(T, k, j);
        double* y = CE(T, k + 1, j);
        double xr = x[0], xi = x[1], yr = y[0], yi = y[1], a, b;
        c_mul(snr, sni, yr, yi, &a, &b);
        x[0] = cs * xr + a; x[1] = cs * xi + b;
        c_mul(snr, -sni, xr, xi, &a, &b);
        y[0] = cs * yr - a; y[1] = cs * yi - b;
    }
    /* columns k, k+1 with (cs, conj(sn)): x' = cs x + conj(sn) y, y' = cs y - sn x */
    for (int pass = 0; pass < 2; pass++) {
        double* M = pass == 0 ? T : Q;
        int rows = pass == 0 ? k : n;
        for (int i = 0; i < rows; i++) {
            double* x = CE(M, i, k);
            double* y = CE(M, i, k + 1);
            double xr = x[0], xi = x[1], yr = y[0], yi = y[1], a, b;
            c_mul(snr, -sni, yr, yi, &a, &b);
            x[0] = cs * xr + a; x[1] = cs * xi + b;
            c_mul(snr, sni, xr, xi, &a, &b);
            y[0] = cs * yr - a; y[1] = cs * yi - b;
        }
    }
    CE(T, k, k)[0] = t22r; CE(T, k, k)[1] = t22i;
    CE(T, k + 1, k + 1)[0] = t11r; CE(T, k + 1, k + 1)[1] = t11i;
}

/* Back-transform: R = Q F Q^H. */
static void back_transform(int n, const double* Q, const double* F, double* R, double* W) {
#ifdef USE_LAPACK
    double one[2] = { 1.0, 0.0 }, zero[2] = { 0.0, 0.0 };
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, (void*)one,
                (const void*)Q, n, (const void*)F, n, (void*)zero, (void*)W, n);
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, n, n, n, (void*)one,
                (const void*)W, n, (const void*)Q, n, (void*)zero, (void*)R, n);
#else
    size_t nn = (size_t)n * (size_t)n;
    mm(n, true, Q, F, W);
    double* QH = (double*)malloc(sizeof(double) * 2 * nn);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            CE(QH, i, j)[0] = CE(Q, j, i)[0];
            CE(QH, i, j)[1] = -CE(Q, j, i)[1];
        }
    mm(n, true, W, QH, R);
    free(QH);
#endif
}

/* Bjorck--Hammarling square root of the upper-triangular T into R.  MF_NOSOL
 * when a zero eigenvalue sits in a nontrivial Jordan structure (r_ii + r_jj
 * == 0 against a nonzero numerator). */
static int sqrtm_tri(int n, const double* T, double* R) {
    memset(R, 0, sizeof(double) * 2 * (size_t)n * (size_t)n);
    for (int j = 0; j < n; j++) {
        c_sqrt(CE(T, j, j)[0], CE(T, j, j)[1], &CE(R, j, j)[0], &CE(R, j, j)[1]);
        for (int i = j - 1; i >= 0; i--) {
            double sr = CE(T, i, j)[0], si = CE(T, i, j)[1];
            for (int k = i + 1; k < j; k++) {
                double pr, pi;
                c_mul(CE(R, i, k)[0], CE(R, i, k)[1], CE(R, k, j)[0], CE(R, k, j)[1], &pr, &pi);
                sr -= pr; si -= pi;
            }
            double dr = CE(R, i, i)[0] + CE(R, j, j)[0];
            double di = CE(R, i, i)[1] + CE(R, j, j)[1];
            if (dr == 0.0 && di == 0.0) {
                if (sr == 0.0 && si == 0.0) { CE(R, i, j)[0] = 0.0; CE(R, i, j)[1] = 0.0; continue; }
                return MF_NOSOL;
            }
            c_div(sr, si, dr, di, &CE(R, i, j)[0], &CE(R, i, j)[1]);
        }
    }
    return all_finite(2 * (size_t)n * (size_t)n, R) ? MF_OK : MF_NOSOL;
}

/* Gauss--Legendre nodes and weights on [0, 1], m points (Newton on P_m). */
static void gauss_legendre01(int m, double* x, double* wt) {
    for (int i = 0; i < m; i++) {
        double z = cos(M_PI * (i + 0.75) / (m + 0.5));
        double pp = 1.0;
        for (int it = 0; it < 100; it++) {
            double p1 = 1.0, p2 = 0.0;
            for (int j = 1; j <= m; j++) {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            pp = m * (z * p1 - p2) / (z * z - 1.0);
            double z1 = z;
            z = z1 - p1 / pp;
            if (fabs(z - z1) < 1e-16) break;
        }
        x[i] = 0.5 * (1.0 - z);
        wt[i] = 1.0 / ((1.0 - z * z) * pp * pp);
    }
}

/* log(T) for upper-triangular T by inverse scaling and squaring.  L receives
 * log(T) (upper triangular). */
static int logm_tri(int n, const double* T0, double* L) {
    static const double theta[8] = { 0.0, 1.586970738772063e-5, 2.313807884242979e-3,
                                     1.938179313533253e-2, 6.209171588994762e-2,
                                     1.276404810806775e-1, 2.060962623452836e-1,
                                     2.879093714241194e-1 };
    size_t nn = (size_t)n * (size_t)n, sz = 2 * nn;
    for (int i = 0; i < n; i++)
        if (CE(T0, i, i)[0] == 0.0 && CE(T0, i, i)[1] == 0.0) return MF_ZEROEIG;

    double* T = (double*)malloc(sizeof(double) * sz);
    double* R = (double*)malloc(sizeof(double) * sz);
    double* Y = (double*)malloc(sizeof(double) * sz);
    if (!T || !R || !Y) { free(T); free(R); free(Y); return MF_FAIL; }
    memcpy(T, T0, sizeof(double) * sz);

    int k = 0, p = 0, m = 7, status = MF_OK;
    for (;;) {
        /* tau = ||T - I||_1 */
        for (size_t q = 0; q < sz; q++) R[q] = T[q];
        for (int i = 0; i < n; i++) CE(R, i, i)[0] -= 1.0;
        double tau = norm1(n, true, R);
        if (tau <= theta[7]) {
            p++;
            int j1 = 7, j2 = 7;
            for (int q = 1; q <= 7; q++) if (tau <= theta[q]) { j1 = q; break; }
            for (int q = 1; q <= 7; q++) if (tau / 2.0 <= theta[q]) { j2 = q; break; }
            if (j1 - j2 <= 1 || p == 2) { m = j1; break; }
        }
        if (k > 120) { status = MF_FAIL; break; }
        int st = sqrtm_tri(n, T, Y);
        if (st != MF_OK) { status = st; break; }
        memcpy(T, Y, sizeof(double) * sz);
        k++;
    }
    if (status == MF_OK) {
        /* R = T - I; L = 2^k sum_j w_j R (I + x_j R)^-1 */
        for (size_t q = 0; q < sz; q++) R[q] = T[q];
        for (int i = 0; i < n; i++) CE(R, i, i)[0] -= 1.0;
        double xg[8], wg[8];
        gauss_legendre01(m, xg, wg);
        memset(L, 0, sizeof(double) * sz);
        for (int g = 0; g < m; g++) {
            /* Solve (I + x_g R) Y = R, column by column (upper triangular). */
            for (int j = 0; j < n; j++) {
                for (int i = j; i >= 0; i--) {
                    double sr = CE(R, i, j)[0], si = CE(R, i, j)[1];
                    for (int q = i + 1; q <= j; q++) {
                        double pr, pi;
                        c_mul(xg[g] * CE(R, i, q)[0], xg[g] * CE(R, i, q)[1],
                              CE(Y, q, j)[0], CE(Y, q, j)[1], &pr, &pi);
                        sr -= pr; si -= pi;
                    }
                    c_div(sr, si, 1.0 + xg[g] * CE(R, i, i)[0], xg[g] * CE(R, i, i)[1],
                          &CE(Y, i, j)[0], &CE(Y, i, j)[1]);
                }
                for (int i = j + 1; i < n; i++) { CE(Y, i, j)[0] = 0.0; CE(Y, i, j)[1] = 0.0; }
            }
            for (size_t q = 0; q < sz; q++) L[q] += wg[g] * Y[q];
        }
        for (size_t q = 0; q < sz; q++) L[q] = ldexp(L[q], k);

        /* Recompute the diagonal and superdiagonal from closed forms. */
        for (int i = 0; i < n; i++)
            c_log(CE(T0, i, i)[0], CE(T0, i, i)[1], &CE(L, i, i)[0], &CE(L, i, i)[1]);
        for (int i = 0; i + 1 < n; i++) {
            double l1r = CE(T0, i, i)[0], l1i = CE(T0, i, i)[1];
            double l2r = CE(T0, i + 1, i + 1)[0], l2i = CE(T0, i + 1, i + 1)[1];
            double tr = CE(T0, i, i + 1)[0], ti = CE(T0, i, i + 1)[1];
            double fr, fi;
            if (l1r == l2r && l1i == l2i) {
                c_div(tr, ti, l1r, l1i, &fr, &fi);
            } else {
                double dr = l2r - l1r, di = l2i - l1i;
                double a1 = c_abs(l1r, l1i), a2 = c_abs(l2r, l2i);
                double qr, qi;
                if (c_abs(dr, di) > 0.25 * (a1 + a2)) {
                    double g1r, g1i, g2r, g2i;
                    c_log(l1r, l1i, &g1r, &g1i);
                    c_log(l2r, l2i, &g2r, &g2i);
                    c_div(g2r - g1r, g2i - g1i, dr, di, &qr, &qi);
                } else {
                    /* (2 atanh(z) + 2 Pi i U) / (l2 - l1), z = (l2 - l1)/(l2 + l1) */
                    double zr, zi, ar, ai, g1r, g1i, g2r, g2i;
                    c_div(dr, di, l2r + l1r, l2i + l1i, &zr, &zi);
                    c_atanh(zr, zi, &ar, &ai);
                    c_log(l1r, l1i, &g1r, &g1i);
                    c_log(l2r, l2i, &g2r, &g2i);
                    double u = ceil(((g2i - g1i) - M_PI) / (2.0 * M_PI));
                    c_div(2.0 * ar, 2.0 * ai + 2.0 * M_PI * u, dr, di, &qr, &qi);
                }
                c_mul(tr, ti, qr, qi, &fr, &fi);
            }
            CE(L, i, i + 1)[0] = fr; CE(L, i, i + 1)[1] = fi;
        }
        if (!all_finite(sz, L)) status = MF_FAIL;
    }
    free(T); free(R); free(Y);
    return status;
}

/* ---- Schur--Parlett ----------------------------------------------------- */

/* Infinity norm (largest row sum of moduli) of an L x L complex column-major
 * block. */
static double blk_norm_inf(int L, const double* B) {
    double best = 0.0;
    for (int i = 0; i < L; i++) {
        double r = 0.0;
        for (int j = 0; j < L; j++) {
            size_t q = (size_t)i + (size_t)j * (size_t)L;
            r += fabs(B[2 * q]) + fabs(B[2 * q + 1]);   /* >= |z|, cheap, no hypot */
        }
        if (r > best) best = r;
    }
    return best;
}

#define MF_ATOM_MAXTERMS 500

/* f on one diagonal block T[s0 .. s0+len-1] by a Taylor series about the mean
 * eigenvalue sigma; the result is written into the same block of F.
 *
 * Truncation follows Higham & Davies (SIAM J. Matrix Anal. Appl. 25(2), 2003,
 * Algorithm 4.4 / MATLAB funm_atom): once a term is negligible relative to the
 * partial sum, the tail is bounded by
 *
 *     ||P_{k+1}|| * mu * max_{0<=j<L} max_z |f^(k+1+j)(z)| / j!,
 *
 * where P_{k+1} = M^{k+1}/(k+1)!, M = T_block - sigma I, mu = ||(I - |N|)^-1 e||
 * with N the strictly upper part of the block, and z ranges over the block's
 * eigenvalues (approximated, as MATLAB does, by sigma and the eigenvalue
 * farthest from it).  The series stops only when that bound is below u ||F||.
 *
 * The earlier rule -- two negligible terms AND at least L terms -- was both
 * slower and wrong: a 100 x 100 random matrix chains ~99 eigenvalues into one
 * block, the forced 99+ terms ran past the derivative cache, and the call
 * failed (MatrixFunction[Cos, m] came back unevaluated). */
static int parlett_atom(int n, const double* T, double* F, int s0, int len,
                        mf_scalar_fn fn, void* ctx) {
    if (len == 1) {
        return fn(ctx, CE(T, s0, s0)[0], CE(T, s0, s0)[1], 0,
                  &CE(F, s0, s0)[0], &CE(F, s0, s0)[1]);
    }
    int L = len;
    size_t ll = (size_t)L * (size_t)L;
    double* M  = (double*)calloc(2 * ll, sizeof(double));
    double* P  = (double*)calloc(2 * ll, sizeof(double));
    double* P2 = (double*)calloc(2 * ll, sizeof(double));
    double* Fb = (double*)calloc(2 * ll, sizeof(double));
    double* y  = (double*)malloc(sizeof(double) * (size_t)L);
    int dcap = MF_ATOM_MAXTERMS + L + 2;
    double* dmax = (double*)malloc(sizeof(double) * (size_t)dcap);   /* derivative bounds */
    if (!M || !P || !P2 || !Fb || !y || !dmax) {
        free(M); free(P); free(P2); free(Fb); free(y); free(dmax);
        return MF_FAIL;
    }
#define BE(B, i, j) ((B) + 2 * ((size_t)(i) + (size_t)(j) * (size_t)L))
    double sr = 0.0, si = 0.0;
    for (int i = 0; i < L; i++) { sr += CE(T, s0 + i, s0 + i)[0]; si += CE(T, s0 + i, s0 + i)[1]; }
    sr /= L; si /= L;
    for (int j = 0; j < L; j++)
        for (int i = 0; i <= j; i++) {
            BE(M, i, j)[0] = CE(T, s0 + i, s0 + j)[0];
            BE(M, i, j)[1] = CE(T, s0 + i, s0 + j)[1];
        }
    for (int i = 0; i < L; i++) { BE(M, i, i)[0] -= sr; BE(M, i, i)[1] -= si; }

    /* The eigenvalue farthest from sigma: the second derivative-bound point. */
    double zfr = sr, zfi = si, far = -1.0;
    for (int i = 0; i < L; i++) {
        double d = c_abs(BE(M, i, i)[0], BE(M, i, i)[1]);
        if (d > far) { far = d; zfr = CE(T, s0 + i, s0 + i)[0]; zfi = CE(T, s0 + i, s0 + i)[1]; }
    }

    /* mu = ||y||_inf with (I - |N|) y = e: unit upper triangular, y >= 1. */
    double mu = 1.0;
    for (int i = L - 1; i >= 0; i--) {
        double acc = 1.0;
        for (int j = i + 1; j < L; j++) acc += c_abs(BE(M, i, j)[0], BE(M, i, j)[1]) * y[j];
        y[i] = acc;
        if (acc > mu) mu = acc;
    }

    int have_d = 0;      /* dmax[0 .. have_d-1] are filled */
    double f0r, f0i;
    int st = fn(ctx, sr, si, 0, &f0r, &f0i);
    bool done = false;
    if (st == MF_OK) {
        for (int i = 0; i < L; i++) { BE(Fb, i, i)[0] = f0r; BE(Fb, i, i)[1] = f0i; }
        memcpy(P, M, sizeof(double) * 2 * ll);               /* P = M^1 / 1! */
        for (int k = 1; k <= MF_ATOM_MAXTERMS; k++) {
            double fnF_old = blk_norm_inf(L, Fb);
            double dr, di;
            st = fn(ctx, sr, si, k, &dr, &di);
            if (st != MF_OK) break;
            for (size_t q = 0; q < ll; q++) {
                double ar, ai;
                c_mul(dr, di, P[2 * q], P[2 * q + 1], &ar, &ai);
                P2[2 * q] = ar; P2[2 * q + 1] = ai;
                Fb[2 * q] += ar; Fb[2 * q + 1] += ai;
            }
            double tn = blk_norm_inf(L, P2);
            /* P <- P M / (k + 1) */
            {
#ifdef USE_LAPACK
                double al[2] = { 1.0 / (double)(k + 1), 0.0 }, zero[2] = { 0.0, 0.0 };
                cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, L, L, L,
                            (void*)al, (const void*)P, L, (const void*)M, L,
                            (void*)zero, (void*)P2, L);
#else
                mm(L, true, P, M, P2);
                for (size_t q = 0; q < 2 * ll; q++) P2[q] /= (double)(k + 1);
#endif
                double* t = P; P = P2; P2 = t;
            }
            double pn = blk_norm_inf(L, P);
            if (pn == 0.0) { done = true; break; }             /* nilpotent: exact */
            if (tn > MF_U * (MF_U + fnF_old)) continue;
            /* Negligible term: bound the tail. */
            int need = k + 1 + L;                              /* orders up to k+L */
            if (need > dcap) need = dcap;
            bool ok = true;
            for (int j = have_d; j < need; j++) {
                double ar, ai, br, bi;
                int s1 = fn(ctx, sr, si, j, &ar, &ai);
                int s2 = (far > 0.0) ? fn(ctx, zfr, zfi, j, &br, &bi) : s1;
                if (s1 != MF_OK || s2 != MF_OK) { ok = false; break; }
                if (far <= 0.0) { br = ar; bi = ai; }
                double a = c_abs(ar, ai), b = c_abs(br, bi);
                dmax[j] = a > b ? a : b;
                have_d = j + 1;
            }
            if (!ok) continue;               /* cannot bound: keep summing */
            double omega = 0.0, jf = 1.0;
            for (int j = 0; j < L && k + 1 + j < have_d; j++) {
                if (j > 0) jf *= (double)j;
                double w = dmax[k + 1 + j] / jf;
                if (w > omega) omega = w;
            }
            double trunc = pn * mu * omega;
            if (isfinite(trunc) && trunc <= MF_U * blk_norm_inf(L, Fb)) { done = true; break; }
        }
        if (st == MF_OK && !done) st = MF_FAIL;
    }
    if (st == MF_OK) {
        for (int j = 0; j < L; j++)
            for (int i = 0; i <= j; i++) {
                CE(F, s0 + i, s0 + j)[0] = BE(Fb, i, j)[0];
                CE(F, s0 + i, s0 + j)[1] = BE(Fb, i, j)[1];
            }
    }
#undef BE
    free(M); free(P); free(P2); free(Fb); free(y); free(dmax);
    return st;
}

/* f(T) by Schur--Parlett; T and Q are reordered in place so that each
 * eigenvalue cluster is contiguous. */
static int schur_parlett(int n, double* T, double* Q, double* F,
                         mf_scalar_fn fn, void* ctx) {
    const double delta = 0.1;
    int* lab = (int*)calloc((size_t)n, sizeof(int));
    int* start = (int*)malloc(sizeof(int) * (size_t)(n + 1));
    if (!lab || !start) { free(lab); free(start); return MF_FAIL; }

    /* Cluster: eigenvalues within delta are chained into one block. */
    int next = 0;
    for (int i = 0; i < n; i++) {
        if (lab[i] == 0) lab[i] = ++next;
        for (int j = i + 1; j < n; j++) {
            if (lab[j] == lab[i]) continue;
            double d = c_abs(CE(T, i, i)[0] - CE(T, j, j)[0], CE(T, i, i)[1] - CE(T, j, j)[1]);
            if (d > delta) continue;
            if (lab[j] == 0) { lab[j] = lab[i]; continue; }
            int lo = lab[i] < lab[j] ? lab[i] : lab[j];
            int hi = lab[i] < lab[j] ? lab[j] : lab[i];
            for (int q = 0; q < n; q++) if (lab[q] == hi) lab[q] = lo;
        }
    }
    /* Reorder: bubble the labels into ascending order with adjacent swaps. */
    for (int pass = 0; pass < n; pass++) {
        bool any = false;
        for (int k = 0; k + 1 < n; k++)
            if (lab[k] > lab[k + 1]) {
                schur_swap(n, T, Q, k);
                int t = lab[k]; lab[k] = lab[k + 1]; lab[k + 1] = t;
                any = true;
            }
        if (!any) break;
    }
    int nb = 0;
    for (int i = 0; i < n; i++) if (i == 0 || lab[i] != lab[i - 1]) start[nb++] = i;
    start[nb] = n;

    memset(F, 0, sizeof(double) * 2 * (size_t)n * (size_t)n);
    int st = MF_OK;
    for (int bj = 0; bj < nb && st == MF_OK; bj++) {
        int sj = start[bj], ej = start[bj + 1];
        st = parlett_atom(n, T, F, sj, ej - sj, fn, ctx);
        if (st != MF_OK) break;
        for (int bi = bj - 1; bi >= 0; bi--) {
            int si_ = start[bi], ei = start[bi + 1];
            /* T_ii X - X T_jj = C, column by column of block j */
            for (int c = sj; c < ej; c++) {
                for (int r = ei - 1; r >= si_; r--) {
                    /* C[r][c] = sum_{q = si_}^{sj-1} F[r][q] T[q][c]
                     *         - sum_{q = ei}^{ej-1}  T[r][q] F[q][c] */
                    double cr = 0.0, ci = 0.0, pr, pi;
                    for (int q = r; q < sj; q++) {
                        c_mul(CE(F, r, q)[0], CE(F, r, q)[1], CE(T, q, c)[0], CE(T, q, c)[1], &pr, &pi);
                        cr += pr; ci += pi;
                    }
                    for (int q = ei; q <= c; q++) {
                        c_mul(CE(T, r, q)[0], CE(T, r, q)[1], CE(F, q, c)[0], CE(F, q, c)[1], &pr, &pi);
                        cr -= pr; ci -= pi;
                    }
                    /* + sum_{q = sj}^{c-1} X[r][q] T[q][c]  (X T_jj, earlier columns) */
                    for (int q = sj; q < c; q++) {
                        c_mul(CE(F, r, q)[0], CE(F, r, q)[1], CE(T, q, c)[0], CE(T, q, c)[1], &pr, &pi);
                        cr += pr; ci += pi;
                    }
                    /* - sum_{q = r+1}^{ei-1} T[r][q] X[q][c]  (T_ii X, later rows) */
                    for (int q = r + 1; q < ei; q++) {
                        c_mul(CE(T, r, q)[0], CE(T, r, q)[1], CE(F, q, c)[0], CE(F, q, c)[1], &pr, &pi);
                        cr -= pr; ci -= pi;
                    }
                    double dr = CE(T, r, r)[0] - CE(T, c, c)[0];
                    double di = CE(T, r, r)[1] - CE(T, c, c)[1];
                    c_div(cr, ci, dr, di, &CE(F, r, c)[0], &CE(F, r, c)[1]);
                }
            }
        }
    }
    free(lab); free(start);
    if (st == MF_OK && !all_finite(2 * (size_t)n * (size_t)n, F)) st = MF_FAIL;
    return st;
}

/* Is f(A) of a real A real?  Every real eigenvalue must give a real value and
 * every complex one a value whose conjugate is its partner's. */
static bool real_consistent(int n, const double* T, const double* F) {
    for (int i = 0; i < n; i++) {
        double lr = CE(T, i, i)[0], li = CE(T, i, i)[1];
        double fr = CE(F, i, i)[0], fi = CE(F, i, i)[1];
        double fa = c_abs(fr, fi);
        if (li == 0.0) {
            if (fabs(fi) > 64.0 * MF_U * fa) return false;
            continue;
        }
        bool found = false;
        double la = c_abs(lr, li);
        for (int j = 0; j < n && !found; j++) {
            if (j == i) continue;
            double mr = CE(T, j, j)[0], mi = CE(T, j, j)[1];
            if (c_abs(mr - lr, mi + li) > 1e-8 * la) continue;
            double gr = CE(F, j, j)[0], gi = CE(F, j, j)[1];
            if (c_abs(gr - fr, gi + fi) <= 1e-6 * (fa + 1e-300)) found = true;
        }
        if (!found) return false;
    }
    return true;
}

/* Is the Schur form T diagonal up to rounding -- a normal matrix (every
 * Hermitian, real symmetric, unitary or skew matrix)?  Then f(A) = Q f(D) Q^H
 * exactly and Schur--Parlett has nothing to do: its clustering made one atom
 * of a Hermitian matrix whose eigenvalues sit within 0.1 of each other and
 * summed a Taylor series over the whole block, 13x slower than the reference
 * for MatrixFunction[Sin, s] of a 200 x 200 SPD s. */
static bool diag_schur(int n, const double* T) {
    double off = 0.0, dmax = 0.0;
    for (int j = 0; j < n; j++) {
        double d = c_abs(CE(T, j, j)[0], CE(T, j, j)[1]);
        if (d > dmax) dmax = d;
        for (int i = 0; i < j; i++) {
            double a = c_abs(CE(T, i, j)[0], CE(T, i, j)[1]);
            if (a > off) off = a;
        }
    }
    return off <= 64.0 * (double)n * MF_U * (dmax > 0.0 ? dmax : 1.0);
}

/* f(T) of a diagonal T: f on each eigenvalue. */
static int diag_fun(int n, const double* T, double* F, mf_scalar_fn fn, void* ctx) {
    memset(F, 0, sizeof(double) * 2 * (size_t)n * (size_t)n);
    for (int i = 0; i < n; i++) {
        int st = fn(ctx, CE(T, i, i)[0], CE(T, i, i)[1], 0, &CE(F, i, i)[0], &CE(F, i, i)[1]);
        if (st != MF_OK) return st;
    }
    return MF_OK;
}

int mf_schur_fun(const double* A, int n, bool cx, int kind,
                 mf_scalar_fn fn, void* ctx, double* F, bool* is_real) {
    *is_real = false;
    size_t nn = (size_t)n * (size_t)n, sz = 2 * nn;
    if (n < 1 || !all_finite((cx ? 2 : 1) * nn, A)) return MF_FAIL;

    if (n == 1) {                       /* the scalar function itself */
        double zr = A[0], zi = cx ? A[1] : 0.0;
        int st;
        if (kind == MF_KIND_LOG)       st = mf_fn_log(NULL, zr, zi, 0, &F[0], &F[1]);
        else if (kind == MF_KIND_SQRT) { c_sqrt(zr, zi, &F[0], &F[1]); st = MF_OK; }
        else                           st = fn(ctx, zr, zi, 0, &F[0], &F[1]);
        if (st != MF_OK) return st;
        if (!cx && F[1] == 0.0) *is_real = true;
        return MF_OK;
    }

    double* T  = (double*)malloc(sizeof(double) * sz);
    double* Q  = (double*)malloc(sizeof(double) * sz);
    double* Fs = (double*)malloc(sizeof(double) * sz);
    double* W  = (double*)malloc(sizeof(double) * sz);
    if (!T || !Q || !Fs || !W) { free(T); free(Q); free(Fs); free(W); return MF_FAIL; }
    int st = schur_complex(A, n, cx, T, Q);
    if (st == MF_OK) {
        if (kind == MF_KIND_LOG)       st = logm_tri(n, T, Fs);
        else if (kind == MF_KIND_SQRT) st = sqrtm_tri(n, T, Fs);
        else if (diag_schur(n, T))     st = diag_fun(n, T, Fs, fn, ctx);
        else                           st = schur_parlett(n, T, Q, Fs, fn, ctx);
    }
    if (st == MF_OK) {
        back_transform(n, Q, Fs, F, W);
        if (!cx && real_consistent(n, T, Fs)) {
            for (size_t k = 0; k < nn; k++) F[2 * k + 1] = 0.0;
            *is_real = true;
        }
        if (!all_finite(sz, F)) st = MF_FAIL;
    }
    free(T); free(Q); free(Fs); free(W);
    return st;
}

/* Mathilda -- restarted Krylov eigensolvers on a matrix-vector callback:
 * a native Krylov-Schur / thick-restart Lanczos, and arpack-ng.
 *
 * See krylov_schur.h for the interface and the algorithm summary. Layout
 * conventions here: the Krylov basis V holds m + 1 vectors of nl = w n
 * doubles each (w = 2 for complex, interleaved), one after another, which is
 * a column-major n x (m + 1) matrix for the real case. The projected matrix
 * H is (m + 1) x m, column-major with leading dimension m + 1, entries w
 * doubles wide. Eigenvectors Y of the m x m Rayleigh quotient are always
 * kept complex (2 doubles per entry, column-major m x m) -- the small
 * problem costs nothing next to the matvecs, and one representation serves
 * the real-symmetric, real-general and complex paths.
 */
#include "krylov_schur.h"
#include "lapack.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_ARPACK
#include <arpack/arpack.h>
#endif

static KsBackend g_backend = KS_BACKEND_AUTO;

KsBackend ks_set_backend(KsBackend b) {
    KsBackend o = g_backend;
    g_backend = b;
    return o;
}

bool ks_arpack_active(void) {
#ifdef USE_ARPACK
    if (g_backend == KS_BACKEND_NATIVE) return false;
    if (g_backend == KS_BACKEND_ARPACK) return true;
    const char* e = getenv("MATHILDA_NO_ARPACK");
    return !(e && *e && strcmp(e, "0") != 0);
#else
    return false;
#endif
}

void ks_result_free(KsResult* r) {
    if (!r) return;
    free(r->re);
    free(r->im);
    free(r->vec);
    memset(r, 0, sizeof *r);
}

/* ---- the deterministic start ---------------------------------------------- */

/* A fixed pseudo-random vector in (-1, 1): a reproducible start that has a
 * component along every eigenvector with probability one (a ramp or a
 * constant vector is orthogonal to whole symmetry classes of eigenvectors
 * of a Toeplitz or Laplacian matrix). */
static void fill_random(double* x, int64_t len, uint64_t* state) {
    for (int64_t i = 0; i < len; i++) {
        *state = *state * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
        x[i] = (double)(*state >> 11) * (2.0 / 9007199254740992.0) - 1.0;
    }
}

/* ---- vector kernels --------------------------------------------------------- */

/* (re, im) = x^H y */
static void v_dot(const double* x, const double* y, int64_t n, bool c, double* re, double* im) {
    if (!c) {
        double s = 0.0;
        for (int64_t i = 0; i < n; i++) s += x[i] * y[i];
        *re = s;
        *im = 0.0;
        return;
    }
    double sr = 0.0, si = 0.0;
    for (int64_t i = 0; i < n; i++) {
        double xr = x[2 * i], xi = x[2 * i + 1], yr = y[2 * i], yi = y[2 * i + 1];
        sr += xr * yr + xi * yi;
        si += xr * yi - xi * yr;
    }
    *re = sr;
    *im = si;
}

/* y += (ar + i ai) x */
static void v_axpy(double ar, double ai, const double* x, double* y, int64_t n, bool c) {
    if (!c) {
        for (int64_t i = 0; i < n; i++) y[i] += ar * x[i];
        return;
    }
    for (int64_t i = 0; i < n; i++) {
        double xr = x[2 * i], xi = x[2 * i + 1];
        y[2 * i] += ar * xr - ai * xi;
        y[2 * i + 1] += ar * xi + ai * xr;
    }
}

static double v_norm(const double* x, int64_t len) {
    /* len doubles: the 2-norm of a real or interleaved complex vector */
    double scale = 0.0, ssq = 1.0;
    for (int64_t i = 0; i < len; i++) {
        double a = fabs(x[i]);
        if (a == 0.0) continue;
        if (scale < a) {
            ssq = 1.0 + ssq * (scale / a) * (scale / a);
            scale = a;
        } else {
            ssq += (a / scale) * (a / scale);
        }
    }
    return scale * sqrt(ssq);
}

static void v_scal(double s, double* x, int64_t len) {
    for (int64_t i = 0; i < len; i++) x[i] *= s;
}

static int g_blas = -1;     /* cblas usable (the LAPACK probe passed) */

static bool have_lapack(void) {
    if (g_blas < 0) g_blas = mathilda_lapack_probe() ? 1 : 0;
    return g_blas == 1;
}

/* Orthogonalise w against the first cnt basis vectors (classical
 * Gram-Schmidt, applied twice), accumulating the coefficients V^H w into h
 * (cnt entries, w-wide). tmp has room for cnt w-wide entries. */
static void orth(const double* V, int64_t nl, int cnt, double* wv, double* h, double* tmp,
                 int64_t n, bool c) {
    int w = c ? 2 : 1;
    memset(h, 0, sizeof(double) * (size_t)w * (size_t)cnt);
    for (int pass = 0; pass < 2; pass++) {
#ifdef USE_LAPACK
        if (!c && have_lapack() && n < INT_MAX) {
            cblas_dgemv(CblasColMajor, CblasTrans, (int)n, cnt, 1.0, V, (int)n, wv, 1, 0.0, tmp, 1);
            cblas_dgemv(CblasColMajor, CblasNoTrans, (int)n, cnt, -1.0, V, (int)n, tmp, 1, 1.0, wv, 1);
            for (int i = 0; i < cnt; i++) h[i] += tmp[i];
            continue;
        }
#endif
        for (int i = 0; i < cnt; i++) v_dot(V + (size_t)i * (size_t)nl, wv, n, c, &tmp[w * i], c ? &tmp[w * i + 1] : &(double){0.0});
        for (int i = 0; i < cnt; i++)
            v_axpy(-tmp[w * i], c ? -tmp[w * i + 1] : 0.0, V + (size_t)i * (size_t)nl, wv, n, c);
        for (int i = 0; i < w * cnt; i++) h[i] += tmp[i];
    }
}

/* X (n x p, w-wide, column-major, ld n) = V (n x m) * Q (m x p, w-wide,
 * column-major, ld m). */
static void basis_times(const double* V, int64_t n, int m, const double* Q, int p, bool c, double* X) {
    int w = c ? 2 : 1;
    int64_t nl = (int64_t)w * n;
#ifdef USE_LAPACK
    if (!c && have_lapack() && n < INT_MAX) {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (int)n, p, m, 1.0, V, (int)n, Q, m,
                    0.0, X, (int)n);
        return;
    }
#endif
    memset(X, 0, sizeof(double) * (size_t)nl * (size_t)p);
    for (int j = 0; j < p; j++)
        for (int i = 0; i < m; i++) {
            const double* q = Q + (size_t)w * ((size_t)i + (size_t)j * (size_t)m);
            if (q[0] == 0.0 && (!c || q[1] == 0.0)) continue;
            v_axpy(q[0], c ? q[1] : 0.0, V + (size_t)i * (size_t)nl, X + (size_t)j * (size_t)nl, n, c);
        }
}

/* ---- the projected eigenproblem --------------------------------------------- */

/* Eigen-decomposition of the m x m leading block of H (ld ldh, w-wide):
 * theta = (re, im), complex eigenvectors Y (m x m, column-major, unit
 * 2-norm). A Hermitian H is symmetrised first and solved with dsyev / zheev
 * (orthonormal Y). Returns false when LAPACK is missing or fails. */
static bool eig_small(const double* H, int ldh, int m, bool c, bool herm,
                      double* re, double* im, double* Y) {
    if (!have_lapack()) return false;
    int w = c ? 2 : 1;
    double* A = malloc(sizeof(double) * (size_t)w * (size_t)m * (size_t)m);
    if (!A) return false;
    for (int j = 0; j < m; j++)
        for (int i = 0; i < m; i++) {
            const double* s = H + (size_t)w * ((size_t)i + (size_t)j * (size_t)ldh);
            double* d = A + (size_t)w * ((size_t)i + (size_t)j * (size_t)m);
            if (herm) {
                const double* t = H + (size_t)w * ((size_t)j + (size_t)i * (size_t)ldh);
                d[0] = 0.5 * (s[0] + t[0]);
                if (c) d[1] = 0.5 * (s[1] - t[1]);
            } else {
                d[0] = s[0];
                if (c) d[1] = s[1];
            }
        }
    bool ok = false;
    if (herm) {
        int info = c ? mat_lapack_zheev(m, A, m, re) : mat_lapack_dsyev(m, A, m, re);
        if (info == 0) {
            for (int j = 0; j < m; j++) im[j] = 0.0;
            for (size_t q = 0; q < (size_t)m * (size_t)m; q++) {
                Y[2 * q] = A[w * q];
                Y[2 * q + 1] = c ? A[w * q + 1] : 0.0;
            }
            ok = true;
        }
    } else if (!c) {
        double* VR = malloc(sizeof(double) * (size_t)m * (size_t)m);
        if (VR && mat_lapack_dgeev(m, A, m, re, im, VR, m) == 0) {
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < m; i++) {
                    double vr, vi;
                    if (im[j] == 0.0) { vr = VR[i + (size_t)j * m]; vi = 0.0; }
                    else if (im[j] > 0.0) { vr = VR[i + (size_t)j * m]; vi = VR[i + (size_t)(j + 1) * m]; }
                    else { vr = VR[i + (size_t)(j - 1) * m]; vi = -VR[i + (size_t)j * m]; }
                    Y[2 * (i + (size_t)j * m)] = vr;
                    Y[2 * (i + (size_t)j * m) + 1] = vi;
                }
            }
            ok = true;
        }
        free(VR);
    } else {
        double* VR = malloc(sizeof(double) * 2 * (size_t)m * (size_t)m);
        double* wv = malloc(sizeof(double) * 2 * (size_t)m);
        if (VR && wv && mat_lapack_zgeev(m, A, m, wv, VR, m) == 0) {
            for (int j = 0; j < m; j++) { re[j] = wv[2 * j]; im[j] = wv[2 * j + 1]; }
            memcpy(Y, VR, sizeof(double) * 2 * (size_t)m * (size_t)m);
            ok = true;
        }
        free(VR);
        free(wv);
    }
    free(A);
    if (!ok) return false;
    /* unit 2-norm columns (LAPACK already gives them; a conjugate-pair
     * column assembled from two real ones needs it) */
    for (int j = 0; j < m; j++) {
        double nr = v_norm(Y + (size_t)2 * (size_t)j * (size_t)m, 2 * (int64_t)m);
        if (nr > 0.0) v_scal(1.0 / nr, Y + (size_t)2 * (size_t)j * (size_t)m, 2 * (int64_t)m);
    }
    return true;
}

/* ---- selection ---------------------------------------------------------------- */

static double goodness(KsWhich w, double re, double im) {
    switch (w) {
        case KS_SM: return -hypot(re, im);
        case KS_LR: return re;
        case KS_SR: return -re;
        case KS_LI: return fabs(im);
        case KS_SI: return -fabs(im);
        default:    return hypot(re, im);
    }
}

/* ord = 0..m-1, best first under `which`. Ties: positive imaginary part
 * first (a conjugate pair stays adjacent, +im leading as in the dense
 * convention), then the larger real part, then the index. KS_BE alternates
 * between the top and the bottom of the real spectrum, top first. */
static void order_ritz(KsWhich which, const double* re, const double* im, int m, int* ord) {
    for (int i = 0; i < m; i++) ord[i] = i;
    if (which == KS_BE) {
        /* ascending by real part, then interleave from both ends */
        int* asc = malloc(sizeof(int) * (size_t)(m > 0 ? m : 1));
        if (!asc) return;
        for (int i = 0; i < m; i++) asc[i] = i;
        for (int i = 1; i < m; i++) {
            int cur = asc[i], j = i;
            while (j > 0 && re[asc[j - 1]] > re[cur]) { asc[j] = asc[j - 1]; j--; }
            asc[j] = cur;
        }
        int lo = 0, hi = m - 1, q = 0;
        while (lo <= hi) {
            ord[q++] = asc[hi--];
            if (lo <= hi) ord[q++] = asc[lo++];
        }
        free(asc);
        return;
    }
    for (int i = 1; i < m; i++) {
        int cur = ord[i], j = i;
        double gc = goodness(which, re[cur], im[cur]);
        while (j > 0) {
            int prev = ord[j - 1];
            double gp = goodness(which, re[prev], im[prev]);
            bool before = gc > gp ||
                (gc == gp && (im[cur] > im[prev] ||
                              (im[cur] == im[prev] && (re[cur] > re[prev] ||
                                                       (re[cur] == re[prev] && cur < prev)))));
            if (!before) break;
            ord[j] = prev;
            j--;
        }
        ord[j] = cur;
    }
}

/* ---- applying OP to a complex vector when OP is real --------------------- */

static bool apply_op(const KsProblem* p, const double* x, bool xc, double* y, double* sr, double* si) {
    if (p->cplx || !xc) return p->op(p->ctx, x, y);
    /* real OP, complex x: OP (xr + i xi) = OP xr + i OP xi */
    int64_t n = p->n;
    for (int64_t i = 0; i < n; i++) sr[i] = x[2 * i];
    if (!p->op(p->ctx, sr, si)) return false;
    for (int64_t i = 0; i < n; i++) y[2 * i] = si[i];
    for (int64_t i = 0; i < n; i++) sr[i] = x[2 * i + 1];
    if (!p->op(p->ctx, sr, si)) return false;
    for (int64_t i = 0; i < n; i++) y[2 * i + 1] = si[i];
    return true;
}

/* ---- native Krylov-Schur ---------------------------------------------------- */

#define HX(H, i, j, ldh, w) ((H) + (size_t)(w) * ((size_t)(i) + (size_t)(j) * (size_t)(ldh)))

/* Ritz vectors x_q = V_m y_{ord[q]} for q < cnt, unit 2-norm, into X (each
 * xw doubles: n for a real Hermitian problem, else 2n). */
static void ritz_vectors(const double* V, int64_t n, int m, bool c, const double* Y, const int* ord,
                         int cnt, bool xcplx, double* X) {
    int w = c ? 2 : 1;
    int64_t nl = (int64_t)w * n;
    int64_t xl = xcplx ? 2 * n : n;
    memset(X, 0, sizeof(double) * (size_t)xl * (size_t)cnt);
    for (int q = 0; q < cnt; q++) {
        const double* y = Y + (size_t)2 * (size_t)ord[q] * (size_t)m;
        double* x = X + (size_t)q * (size_t)xl;
        for (int j = 0; j < m; j++) {
            double yr = y[2 * j], yi = y[2 * j + 1];
            if (yr == 0.0 && yi == 0.0) continue;
            const double* v = V + (size_t)j * (size_t)nl;
            if (c) {
                v_axpy(yr, yi, v, x, n, true);
            } else if (xcplx) {
                for (int64_t i = 0; i < n; i++) { x[2 * i] += yr * v[i]; x[2 * i + 1] += yi * v[i]; }
            } else {
                for (int64_t i = 0; i < n; i++) x[i] += yr * v[i];
            }
        }
        double nr = v_norm(x, xl);
        if (nr > 0.0) v_scal(1.0 / nr, x, xl);
    }
}

static KsStatus native_solve(const KsProblem* p, int m, int64_t maxit, double tol, KsResult* r) {
    if (!have_lapack()) return KS_FAIL;
    const int64_t n = p->n;
    const bool c = p->cplx, herm = p->hermitian;
    const int w = c ? 2 : 1;
    const int64_t nl = (int64_t)w * n;
    const int nev = p->nev;
    const int ldh = m + 1;
    /* Ritz vectors come out complex unless OP is real and Hermitian */
    const bool xcplx = c || !herm;
    const int64_t xl = xcplx ? 2 * n : n;
    const double eps = DBL_EPSILON;
    const double eps23 = pow(eps, 2.0 / 3.0);

    double* V = calloc((size_t)nl * (size_t)(m + 1), sizeof(double));
    double* H = calloc((size_t)w * (size_t)ldh * (size_t)m, sizeof(double));
    double* Y = malloc(sizeof(double) * 2 * (size_t)m * (size_t)m);
    double* tre = malloc(sizeof(double) * (size_t)m);
    double* tim = malloc(sizeof(double) * (size_t)m);
    int* ord = malloc(sizeof(int) * (size_t)m);
    double* h = malloc(sizeof(double) * (size_t)w * (size_t)(m + 1));
    double* tmp = malloc(sizeof(double) * (size_t)w * (size_t)(m + 1));
    double* Q = malloc(sizeof(double) * (size_t)w * (size_t)m * (size_t)m);
    double* T = malloc(sizeof(double) * (size_t)w * (size_t)m * (size_t)m);
    double* Hn = malloc(sizeof(double) * (size_t)w * (size_t)m * (size_t)m);
    double* Vt = malloc(sizeof(double) * (size_t)nl * (size_t)m);
    double* X = malloc(sizeof(double) * (size_t)xl * (size_t)nev);
    double* ox = malloc(sizeof(double) * (size_t)xl);
    double* s1 = malloc(sizeof(double) * (size_t)n);
    double* s2 = malloc(sizeof(double) * (size_t)n);
    KsStatus st = KS_NOMEM;
    if (!V || !H || !Y || !tre || !tim || !ord || !h || !tmp || !Q || !T || !Hn || !Vt || !X ||
        !ox || !s1 || !s2)
        goto done;

    uint64_t seed = UINT64_C(0x9E3779B97F4A7C15);
    if (p->v0) memcpy(V, p->v0, sizeof(double) * (size_t)nl);
    else fill_random(V, nl, &seed);
    {
        double nr = v_norm(V, nl);
        if (!(nr > 0.0) || !isfinite(nr)) { fill_random(V, nl, &seed); nr = v_norm(V, nl); }
        v_scal(1.0 / nr, V, nl);
    }

    st = KS_MAXIT;
    int k = 0;
    int64_t it;
    for (it = 0; it < maxit; it++) {
        /* ---- expand the decomposition from k to m vectors ---- */
        for (int j = k; j < m; j++) {
            double* vj = V + (size_t)j * (size_t)nl;
            double* wv = V + (size_t)(j + 1) * (size_t)nl;
            if (!p->op(p->ctx, vj, wv)) { st = KS_FAIL; goto done; }
            double wn0 = v_norm(wv, nl);
            if (!isfinite(wn0)) { st = KS_FAIL; goto done; }
            orth(V, nl, j + 1, wv, h, tmp, n, c);
            for (int i = 0; i <= j; i++) {
                double* hij = HX(H, i, j, ldh, w);
                hij[0] = h[w * i];
                if (c) hij[1] = h[w * i + 1];
            }
            double beta = v_norm(wv, nl);
            double* hb = HX(H, j + 1, j, ldh, w);
            if (beta <= 1e-12 * wn0 || beta == 0.0) {
                /* an invariant subspace: the Ritz values in it are exact;
                 * continue with a fresh direction orthogonal to it */
                hb[0] = 0.0;
                if (c) hb[1] = 0.0;
                if (j + 1 < n) {
                    fill_random(wv, nl, &seed);
                    orth(V, nl, j + 1, wv, h, tmp, n, c);
                    double nr = v_norm(wv, nl);
                    if (nr > 1e-8) v_scal(1.0 / nr, wv, nl);
                    else memset(wv, 0, sizeof(double) * (size_t)nl);
                } else {
                    memset(wv, 0, sizeof(double) * (size_t)nl);
                }
            } else {
                v_scal(1.0 / beta, wv, nl);
                hb[0] = beta;
                if (c) hb[1] = 0.0;
            }
        }

        /* ---- Ritz pairs of the projected matrix ---- */
        if (!eig_small(H, ldh, m, c, herm, tre, tim, Y)) { st = KS_FAIL; goto done; }
        order_ritz(p->which, tre, tim, m, ord);
        double betaL = HX(H, m, m - 1, ldh, w)[0];
        double tmax = 0.0;
        for (int i = 0; i < m; i++) tmax = fmax(tmax, hypot(tre[i], tim[i]));
        int conv = 0;
        for (int q = 0; q < nev; q++) {
            int i = ord[q];
            const double* y = Y + (size_t)2 * (size_t)i * (size_t)m;
            double est = betaL * hypot(y[2 * (m - 1)], y[2 * (m - 1) + 1]);
            if (est <= tol * fmax(hypot(tre[i], tim[i]), eps23 * tmax)) conv++;
        }
        if (conv >= nev) {
            /* confirm by the explicit residual ||OP x - theta x|| */
            ritz_vectors(V, n, m, c, Y, ord, nev, xcplx, X);
            bool good = true;
            for (int q = 0; q < nev && good; q++) {
                int i = ord[q];
                const double* x = X + (size_t)q * (size_t)xl;
                if (!apply_op(p, x, xcplx, ox, s1, s2)) { st = KS_FAIL; goto done; }
                if (xcplx) v_axpy(-tre[i], -tim[i], x, ox, n, true);
                else v_axpy(-tre[i], 0.0, x, ox, n, false);
                double res = v_norm(ox, xl);
                double thr = fmax(100.0 * tol * fmax(hypot(tre[i], tim[i]), eps23 * tmax),
                                  1e3 * eps * tmax);
                if (!(res <= thr)) good = false;
            }
            if (good) { st = KS_OK; it++; break; }
        }
        if (it + 1 >= maxit) { it++; break; }
        if (betaL == 0.0 && m >= n) { st = KS_OK; it++; break; }   /* the whole space */

        /* ---- thick restart onto the wanted invariant subspace ---- */
        int keep = nev + (m - nev) / 2;
        if (keep > m - 1) keep = m - 1;
        if (keep < 1) keep = 1;
        /* Q: the basis columns, real for a real OP */
        int cols = 0;
        if (herm) {
            for (int q = 0; q < keep; q++) {
                const double* y = Y + (size_t)2 * (size_t)ord[q] * (size_t)m;
                double* qc = Q + (size_t)w * (size_t)cols * (size_t)m;
                for (int i = 0; i < m; i++) {
                    qc[w * i] = y[2 * i];
                    if (c) qc[w * i + 1] = y[2 * i + 1];
                }
                cols++;
            }
        } else {
            for (int q = 0; q < keep && cols < m - 1; q++) {
                int i = ord[q];
                const double* y = Y + (size_t)2 * (size_t)i * (size_t)m;
                int add = (!c && tim[i] != 0.0) ? 2 : 1;
                if (!c && tim[i] != 0.0) {
                    /* the conjugate partner shares the span: skip it */
                    bool seen = false;
                    for (int q2 = 0; q2 < q; q2++) {
                        int i2 = ord[q2];
                        if (tim[i2] == -tim[i] && tre[i2] == tre[i]) { seen = true; break; }
                    }
                    if (seen) continue;
                    if (cols + add > m - 1) break;
                }
                for (int part = 0; part < add; part++) {
                    double* qc = Q + (size_t)w * (size_t)cols * (size_t)m;
                    for (int r2 = 0; r2 < m; r2++) {
                        if (c) { qc[2 * r2] = y[2 * r2]; qc[2 * r2 + 1] = y[2 * r2 + 1]; }
                        else qc[r2] = y[2 * r2 + part];
                    }
                    /* modified Gram-Schmidt, twice */
                    for (int pass = 0; pass < 2; pass++)
                        for (int prev = 0; prev < cols; prev++) {
                            double dr, di;
                            const double* qp = Q + (size_t)w * (size_t)prev * (size_t)m;
                            v_dot(qp, qc, m, c, &dr, &di);
                            v_axpy(-dr, -di, qp, qc, m, c);
                        }
                    double nr = v_norm(qc, (int64_t)w * m);
                    if (nr > 1e-10) { v_scal(1.0 / nr, qc, (int64_t)w * m); cols++; }
                }
            }
            if (cols < 1) { st = KS_FAIL; goto done; }
        }
        /* H_new = Q^H H_m Q (a Hermitian OP: diag(theta) exactly) */
        if (herm) {
            memset(Hn, 0, sizeof(double) * (size_t)w * (size_t)cols * (size_t)cols);
            for (int q = 0; q < cols; q++) HX(Hn, q, q, cols, w)[0] = tre[ord[q]];
        } else {
            /* T = H_m Q */
            memset(T, 0, sizeof(double) * (size_t)w * (size_t)m * (size_t)cols);
            for (int jj = 0; jj < cols; jj++)
                for (int l = 0; l < m; l++) {
                    const double* ql = Q + (size_t)w * ((size_t)l + (size_t)jj * (size_t)m);
                    if (ql[0] == 0.0 && (!c || ql[1] == 0.0)) continue;
                    for (int i = 0; i < m; i++) {
                        const double* hil = HX(H, i, l, ldh, w);
                        double* t = T + (size_t)w * ((size_t)i + (size_t)jj * (size_t)m);
                        if (c) {
                            t[0] += hil[0] * ql[0] - hil[1] * ql[1];
                            t[1] += hil[0] * ql[1] + hil[1] * ql[0];
                        } else {
                            t[0] += hil[0] * ql[0];
                        }
                    }
                }
            for (int jj = 0; jj < cols; jj++)
                for (int ii = 0; ii < cols; ii++) {
                    double dr, di;
                    v_dot(Q + (size_t)w * (size_t)ii * (size_t)m, T + (size_t)w * (size_t)jj * (size_t)m,
                          m, c, &dr, &di);
                    double* hn = HX(Hn, ii, jj, cols, w);
                    hn[0] = dr;
                    if (c) hn[1] = di;
                }
        }
        /* V_new = V_m Q; the residual direction follows */
        basis_times(V, n, m, Q, cols, c, Vt);
        memcpy(V, Vt, sizeof(double) * (size_t)nl * (size_t)cols);
        memmove(V + (size_t)cols * (size_t)nl, V + (size_t)m * (size_t)nl, sizeof(double) * (size_t)nl);
        /* the new projected matrix and its spike row beta e_m^T Q */
        memset(H, 0, sizeof(double) * (size_t)w * (size_t)ldh * (size_t)m);
        for (int jj = 0; jj < cols; jj++) {
            for (int ii = 0; ii < cols; ii++) {
                const double* s = HX(Hn, ii, jj, cols, w);
                double* d = HX(H, ii, jj, ldh, w);
                d[0] = s[0];
                if (c) d[1] = s[1];
            }
            const double* ql = Q + (size_t)w * ((size_t)(m - 1) + (size_t)jj * (size_t)m);
            double* d = HX(H, cols, jj, ldh, w);
            d[0] = betaL * ql[0];
            if (c) d[1] = betaL * ql[1];
        }
        k = cols;
    }

    /* ---- the answer: the nev wanted Ritz pairs of the last expansion ---- */
    r->re = malloc(sizeof(double) * (size_t)nev);
    r->im = malloc(sizeof(double) * (size_t)nev);
    if (!r->re || !r->im) { ks_result_free(r); st = KS_NOMEM; goto done; }
    for (int q = 0; q < nev; q++) {
        r->re[q] = tre[ord[q]];
        r->im[q] = herm ? 0.0 : tim[ord[q]];
    }
    if (p->want_vectors) {
        ritz_vectors(V, n, m, c, Y, ord, nev, xcplx, X);
        r->vec = X;
        X = NULL;
        r->vec_cplx = xcplx;
    }
    r->nconv = nev;
    r->iters = it;
    r->backend = "KrylovSchur";

done:
    free(V); free(H); free(Y); free(tre); free(tim); free(ord); free(h); free(tmp);
    free(Q); free(T); free(Hn); free(Vt); free(X); free(ox); free(s1); free(s2);
    return st;
}

/* ---- ARPACK ------------------------------------------------------------------- */

#ifdef USE_ARPACK

static const char* arpack_which(KsWhich w, bool sym) {
    switch (w) {
        case KS_SM: return "SM";
        case KS_LR: return sym ? "LA" : "LR";
        case KS_SR: return sym ? "SA" : "SR";
        case KS_LI: return sym ? "LM" : "LI";
        case KS_SI: return sym ? "SM" : "SI";
        case KS_BE: return sym ? "BE" : NULL;
        default:    return "LM";
    }
}

/* Gather nconv pairs (vals re/im, vectors xl doubles each) into r, best
 * first under `which`, keeping nev. Takes ownership of nothing. */
static KsStatus arpack_finish(const KsProblem* p, int nconv, const double* vr, const double* vi,
                              const double* X, bool xcplx, KsResult* r) {
    int nev = p->nev;
    if (nconv < nev) return KS_FAIL;
    int* ord = malloc(sizeof(int) * (size_t)nconv);
    if (!ord) return KS_NOMEM;
    order_ritz(p->which, vr, vi, nconv, ord);
    int64_t xl = xcplx ? 2 * p->n : p->n;
    r->re = malloc(sizeof(double) * (size_t)nev);
    r->im = malloc(sizeof(double) * (size_t)nev);
    if (p->want_vectors) r->vec = malloc(sizeof(double) * (size_t)xl * (size_t)nev);
    if (!r->re || !r->im || (p->want_vectors && !r->vec)) { free(ord); ks_result_free(r); return KS_NOMEM; }
    for (int q = 0; q < nev; q++) {
        r->re[q] = vr[ord[q]];
        r->im[q] = p->hermitian ? 0.0 : vi[ord[q]];
        if (p->want_vectors) {
            double* x = r->vec + (size_t)q * (size_t)xl;
            memcpy(x, X + (size_t)ord[q] * (size_t)xl, sizeof(double) * (size_t)xl);
            double nr = v_norm(x, xl);
            if (nr > 0.0) v_scal(1.0 / nr, x, xl);
        }
    }
    free(ord);
    r->vec_cplx = xcplx;
    r->nconv = nev;
    r->backend = "ARPACK";
    return KS_OK;
}

static KsStatus arpack_solve(const KsProblem* p, int ncv, int64_t maxit, double tol, KsResult* r) {
    const a_int n = (a_int)p->n;
    const a_int nev = p->nev;
    const bool c = p->cplx, sym = !c && p->hermitian;
    const char* which = arpack_which(p->which, sym);   /* znaupd takes the general names */
    if (!which) return KS_FAIL;
    /* the library's own limits */
    if (sym) { if (!(nev < n && ncv > nev && ncv <= n)) return KS_FAIL; }
    else { if (!(nev <= n - 2 && ncv - nev >= 2 && ncv <= n)) return KS_FAIL; }
    const int w = c ? 2 : 1;
    a_int iparam[11] = { 0 }, ipntr[14] = { 0 };
    iparam[0] = 1;                       /* exact shifts */
    iparam[2] = (a_int)maxit;
    iparam[6] = 1;                       /* mode 1: OP x = lambda x */
    a_int lworkl = sym ? ncv * (ncv + 8) : (c ? 3 * ncv * ncv + 5 * ncv : 3 * ncv * ncv + 6 * ncv);
    double* resid = malloc(sizeof(double) * (size_t)w * (size_t)n);
    double* V = malloc(sizeof(double) * (size_t)w * (size_t)n * (size_t)ncv);
    double* workd = malloc(sizeof(double) * (size_t)w * 3 * (size_t)n);
    double* workl = malloc(sizeof(double) * (size_t)w * (size_t)lworkl);
    double* rwork = c ? malloc(sizeof(double) * (size_t)ncv) : NULL;
    KsStatus st = KS_NOMEM;
    double *d = NULL, *di = NULL, *Z = NULL, *workev = NULL, *X = NULL;
    a_int* sel = NULL;
    if (!resid || !V || !workd || !workl || (c && !rwork)) goto out;
    uint64_t seed = UINT64_C(0x9E3779B97F4A7C15);
    if (p->v0) memcpy(resid, p->v0, sizeof(double) * (size_t)w * (size_t)n);
    else fill_random(resid, (int64_t)w * n, &seed);
    a_int ido = 0, info = 1;
    for (;;) {
        if (sym)
            dsaupd_c(&ido, "I", n, which, nev, tol, resid, ncv, V, n, iparam, ipntr, workd, workl, lworkl, &info);
        else if (!c)
            dnaupd_c(&ido, "I", n, which, nev, tol, resid, ncv, V, n, iparam, ipntr, workd, workl, lworkl, &info);
        else
            znaupd_c(&ido, "I", n, which, nev, tol, (a_dcomplex*)resid, ncv, (a_dcomplex*)V, n, iparam, ipntr,
                     (a_dcomplex*)workd, (a_dcomplex*)workl, lworkl, rwork, &info);
        if (ido != -1 && ido != 1) break;
        const double* x = workd + (size_t)w * (size_t)(ipntr[0] - 1);
        double* y = workd + (size_t)w * (size_t)(ipntr[1] - 1);
        if (!p->op(p->ctx, x, y)) { st = KS_FAIL; goto out; }
    }
    if (info != 0) { st = KS_FAIL; goto out; }     /* maxit (1), no shifts (3), errors: native */
    a_int nconv = iparam[4];
    if (nconv < nev) { st = KS_FAIL; goto out; }
    a_int rvec = p->want_vectors ? 1 : 0;
    sel = calloc((size_t)ncv, sizeof(a_int));
    if (!sel) goto out;
    if (sym) {
        d = malloc(sizeof(double) * (size_t)nev);
        di = calloc((size_t)nev, sizeof(double));
        Z = rvec ? malloc(sizeof(double) * (size_t)n * (size_t)nev) : NULL;
        if (!d || !di || (rvec && !Z)) goto out;
        dseupd_c(rvec, "A", sel, d, rvec ? Z : V, n, 0.0, "I", n, which, nev, tol, resid, ncv, V, n,
                 iparam, ipntr, workd, workl, lworkl, &info);
        if (info != 0) { st = KS_FAIL; goto out; }
        st = arpack_finish(p, (int)nconv, d, di, Z, false, r);
    } else if (!c) {
        d = malloc(sizeof(double) * (size_t)(nev + 1));
        di = malloc(sizeof(double) * (size_t)(nev + 1));
        Z = malloc(sizeof(double) * (size_t)n * (size_t)(nev + 1));
        workev = malloc(sizeof(double) * 3 * (size_t)ncv);
        if (!d || !di || !Z || !workev) goto out;
        dneupd_c(rvec, "A", sel, d, di, Z, n, 0.0, 0.0, workev, "I", n, which, nev, tol, resid, ncv, V, n,
                 iparam, ipntr, workd, workl, lworkl, &info);
        if (info != 0) { st = KS_FAIL; goto out; }
        nconv = iparam[4];
        if (nconv > nev + 1) nconv = nev + 1;
        if (rvec) {
            /* unpack dgeev-style conjugate pairs into complex vectors */
            X = malloc(sizeof(double) * 2 * (size_t)n * (size_t)nconv);
            if (!X) goto out;
            for (a_int j = 0; j < nconv; j++) {
                double* x = X + (size_t)2 * (size_t)n * (size_t)j;
                const double* zr;
                const double* zi = NULL;
                double sgn = 1.0;
                if (di[j] == 0.0) zr = Z + (size_t)n * (size_t)j;
                else if (di[j] > 0.0 && j + 1 <= nev) { zr = Z + (size_t)n * (size_t)j; zi = Z + (size_t)n * (size_t)(j + 1); }
                else if (di[j] < 0.0 && j >= 1) { zr = Z + (size_t)n * (size_t)(j - 1); zi = Z + (size_t)n * (size_t)j; sgn = -1.0; }
                else { st = KS_FAIL; goto out; }
                for (a_int i = 0; i < n; i++) { x[2 * i] = zr[i]; x[2 * i + 1] = zi ? sgn * zi[i] : 0.0; }
            }
        }
        st = arpack_finish(p, (int)nconv, d, di, X, true, r);
    } else {
        d = malloc(sizeof(double) * 2 * (size_t)(nev + 1));
        di = malloc(sizeof(double) * (size_t)(nev + 1));
        Z = rvec ? malloc(sizeof(double) * 2 * (size_t)n * (size_t)nev) : NULL;
        workev = malloc(sizeof(double) * 4 * (size_t)ncv);
        if (!d || !di || (rvec && !Z) || !workev) goto out;
        a_dcomplex sigma;
        memset(&sigma, 0, sizeof sigma);
        zneupd_c(rvec, "A", sel, (a_dcomplex*)d, (a_dcomplex*)(rvec ? Z : V), n, sigma, (a_dcomplex*)workev,
                 "I", n, which, nev, tol, (a_dcomplex*)resid, ncv, (a_dcomplex*)V, n, iparam, ipntr,
                 (a_dcomplex*)workd, (a_dcomplex*)workl, lworkl, rwork, &info);
        if (info != 0) { st = KS_FAIL; goto out; }
        if (nconv > nev) nconv = nev;
        double* vre = malloc(sizeof(double) * (size_t)nconv);
        if (!vre) goto out;
        for (a_int j = 0; j < nconv; j++) { vre[j] = d[2 * j]; di[j] = d[2 * j + 1]; }
        st = arpack_finish(p, (int)nconv, vre, di, Z, true, r);
        free(vre);
    }
    if (st == KS_OK) r->iters = iparam[2];
out:
    free(resid); free(V); free(workd); free(workl); free(rwork);
    free(d); free(di); free(Z); free(workev); free(X); free(sel);
    return st;
}
#endif /* USE_ARPACK */

/* ---- entry ------------------------------------------------------------------ */

KsStatus ks_solve(const KsProblem* p, KsResult* r) {
    memset(r, 0, sizeof *r);
    if (!p || !p->op || p->n < 1 || p->nev < 1 || p->nev > p->n) return KS_FAIL;
    if (p->which == KS_BE && !p->hermitian) return KS_FAIL;
    int64_t n = p->n;
    int ncv = p->ncv;
    if (ncv <= 0) {
        int64_t d = 2 * (int64_t)p->nev + 1;
        if (d < 20) d = 20;
        ncv = (int)(d < n ? d : n);
    }
    if (ncv > n) ncv = (int)n;
    if (ncv < p->nev) ncv = p->nev;
    int64_t maxit = p->maxit > 0 ? p->maxit : 1000;
    double tol = p->tol > 0.0 ? p->tol : DBL_EPSILON;
#ifdef USE_ARPACK
    if (ks_arpack_active()) {
        KsStatus st = arpack_solve(p, ncv, maxit, tol, r);
        if (st == KS_OK || st == KS_NOMEM) return st;
        ks_result_free(r);
    }
#endif
    /* the native basis needs one spare column past nev to restart; when
     * ncv == nev < n, widen it */
    if (ncv == p->nev && ncv < n) ncv++;
    return native_solve(p, ncv, maxit, tol, r);
}

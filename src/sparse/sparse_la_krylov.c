/* Mathilda -- sparse linear algebra: Krylov solvers and preconditioners.
 *
 *   CG        Hermitian positive definite A (preconditioned conjugate gradient)
 *   BiCGSTAB  general A (van der Vorst), right-preconditioned  -- the default,
 *             with ILU0, as the reference's Method -> "Krylov" default
 *   GMRES(m)  general A, restarted, right-preconditioned, complex Givens
 *
 * Preconditioners: None, Jacobi, ILU0 (no fill, on A's own pattern), ILUT
 * (Saad's dual-threshold incomplete LU) and ILUTP (ILUT with column
 * pivoting), or a user callback.
 *
 * Real and complex share one code path: vectors are double arrays with
 * w = 1 (real) or w = 2 (interleaved complex) doubles per entry, and every
 * scalar is a Cx {re, im} whose imaginary part stays 0 in the real case.
 * Inner products conjugate their first argument.
 */
#include "sparse_la.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct { double re, im; } Cx;

static inline Cx cx(double r, double i) { Cx z = { r, i }; return z; }
static inline Cx cmul(Cx a, Cx b) { return cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
static inline Cx cdiv(Cx a, Cx b) {
    double d = b.re * b.re + b.im * b.im;
    return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}
static inline Cx cconj(Cx a) { return cx(a.re, -a.im); }
static inline double cabs2(Cx a) { return a.re * a.re + a.im * a.im; }

typedef struct {
    int64_t n;
    int w;          /* 1 real, 2 complex */
} Vsp;

static inline Cx vget(const Vsp* s, const double* v, int64_t i) {
    return s->w == 1 ? cx(v[i], 0.0) : cx(v[2 * i], v[2 * i + 1]);
}
static inline void vset(const Vsp* s, double* v, int64_t i, Cx z) {
    if (s->w == 1) v[i] = z.re;
    else { v[2 * i] = z.re; v[2 * i + 1] = z.im; }
}

/* <a, b> = sum conj(a_i) b_i */
static Cx vdot(const Vsp* s, const double* a, const double* b) {
    if (s->w == 1) {
        double acc = 0.0;
        for (int64_t i = 0; i < s->n; i++) acc += a[i] * b[i];
        return cx(acc, 0.0);
    }
    double re = 0.0, im = 0.0;
    for (int64_t i = 0; i < s->n; i++) {
        double ar = a[2 * i], ai = a[2 * i + 1], br = b[2 * i], bi = b[2 * i + 1];
        re += ar * br + ai * bi;
        im += ar * bi - ai * br;
    }
    return cx(re, im);
}
static double vnorm(const Vsp* s, const double* a) { return sqrt(vdot(s, a, a).re); }
/* y += alpha x */
static void vaxpy(const Vsp* s, Cx alpha, const double* x, double* y) {
    if (s->w == 1) {
        for (int64_t i = 0; i < s->n; i++) y[i] += alpha.re * x[i];
        return;
    }
    for (int64_t i = 0; i < s->n; i++) {
        double xr = x[2 * i], xi = x[2 * i + 1];
        y[2 * i] += alpha.re * xr - alpha.im * xi;
        y[2 * i + 1] += alpha.re * xi + alpha.im * xr;
    }
}
static void vscale(const Vsp* s, Cx alpha, double* x) {
    for (int64_t i = 0; i < s->n; i++) vset(s, x, i, cmul(alpha, vget(s, x, i)));
}
static void vcopy(const Vsp* s, const double* x, double* y) {
    memcpy(y, x, sizeof(double) * (size_t)s->w * (size_t)s->n);
}
static double* vnew(const Vsp* s) {
    return calloc((size_t)s->w * (size_t)(s->n > 0 ? s->n : 1), sizeof(double));
}

/* ======================================================================== */
/* Preconditioners                                                          */
/* ======================================================================== */

/* An incomplete factorisation P A Q ~ L U kept as two CSR matrices: L strict
 * lower (unit diagonal implied), U upper with the diagonal FIRST in each row
 * (ud[i] holds its index). `perm` (ILUTP) maps the factor's column j to A's
 * column perm[j]; NULL for no column permutation. */
typedef struct {
    KryPrec kind;
    Vsp s;
    int64_t *lp, *li, *up, *ui;
    double *lx, *ux;
    int64_t* perm;
    double* dinv;       /* Jacobi */
    KryUserPrec user;
    void* user_ctx;
    double* tmp;
} Prec;

static void prec_free(Prec* p) {
    free(p->lp); free(p->li); free(p->lx);
    free(p->up); free(p->ui); free(p->ux);
    free(p->perm); free(p->dinv); free(p->tmp);
    memset(p, 0, sizeof *p);
}

/* z = M^-1 r */
static bool prec_apply(const Prec* p, const double* r, double* z) {
    const Vsp* s = &p->s;
    switch (p->kind) {
        case KPC_NONE:
            vcopy(s, r, z);
            return true;
        case KPC_USER:
            return p->user(p->user_ctx, r, z);
        case KPC_JACOBI:
            for (int64_t i = 0; i < s->n; i++)
                vset(s, z, i, cmul(vget(s, p->dinv, i), vget(s, r, i)));
            return true;
        default: break;
    }
    /* L y = r, U t = y, z = Q t */
    double* y = p->tmp;
    for (int64_t i = 0; i < s->n; i++) {
        Cx acc = vget(s, r, i);
        for (int64_t k = p->lp[i]; k < p->lp[i + 1]; k++) {
            Cx l = s->w == 1 ? cx(p->lx[k], 0.0) : cx(p->lx[2 * k], p->lx[2 * k + 1]);
            Cx yj = vget(s, y, p->li[k]);
            acc.re -= l.re * yj.re - l.im * yj.im;
            acc.im -= l.re * yj.im + l.im * yj.re;
        }
        vset(s, y, i, acc);
    }
    for (int64_t i = s->n; i-- > 0;) {
        Cx acc = vget(s, y, i);
        int64_t d = p->up[i];
        for (int64_t k = d + 1; k < p->up[i + 1]; k++) {
            Cx u = s->w == 1 ? cx(p->ux[k], 0.0) : cx(p->ux[2 * k], p->ux[2 * k + 1]);
            Cx yj = vget(s, y, p->ui[k]);
            acc.re -= u.re * yj.re - u.im * yj.im;
            acc.im -= u.re * yj.im + u.im * yj.re;
        }
        Cx dg = s->w == 1 ? cx(p->ux[d], 0.0) : cx(p->ux[2 * d], p->ux[2 * d + 1]);
        vset(s, y, i, cdiv(acc, dg));
    }
    if (p->perm) {
        for (int64_t j = 0; j < s->n; j++) vset(s, z, p->perm[j], vget(s, y, j));
    } else {
        vcopy(s, y, z);
    }
    return true;
}

static inline Cx mget(const SpMat* A, int64_t k) {
    return A->cplx ? cx(A->x[2 * k], A->x[2 * k + 1]) : cx(A->x[k], 0.0);
}

static bool prec_jacobi(const SpMat* A, Prec* p) {
    p->dinv = vnew(&p->s);
    if (!p->dinv) return false;
    for (int64_t i = 0; i < A->n; i++) {
        Cx d = cx(0.0, 0.0);
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++)
            if (A->ci[k] == i) d = mget(A, k);
        vset(&p->s, p->dinv, i, cabs2(d) > 0.0 ? cdiv(cx(1.0, 0.0), d) : cx(1.0, 0.0));
    }
    return true;
}

/* Growable CSR row builder for L and U. */
typedef struct { int64_t *p, *i; double* x; int64_t len, cap; int w; } Grow;

static bool grow_push(Grow* g, int64_t col, Cx v) {
    if (g->len == g->cap) {
        int64_t nc = g->cap ? 2 * g->cap : 1024;
        int64_t* ni = realloc(g->i, sizeof(int64_t) * (size_t)nc);
        if (!ni) return false;
        g->i = ni;
        double* nx = realloc(g->x, sizeof(double) * (size_t)g->w * (size_t)nc);
        if (!nx) return false;
        g->x = nx;
        g->cap = nc;
    }
    g->i[g->len] = col;
    g->x[g->w * g->len] = v.re;
    if (g->w == 2) g->x[2 * g->len + 1] = v.im;
    g->len++;
    return true;
}

/* ILU(0): the factors keep exactly A's pattern (IKJ elimination on a copy).
 * Fails on a structurally missing or zero pivot. */
static bool prec_ilu0(const SpMat* A, Prec* p) {
    int64_t n = A->n;
    int w = p->s.w;
    Grow L = { 0, 0, 0, 0, 0, w }, U = { 0, 0, 0, 0, 0, w };
    p->lp = malloc(sizeof(int64_t) * (size_t)(n + 1));
    p->up = malloc(sizeof(int64_t) * (size_t)(n + 1));
    int64_t* where = malloc(sizeof(int64_t) * (size_t)(n > 0 ? n : 1));
    double* row = vnew(&p->s);
    if (!p->lp || !p->up || !where || !row) { free(where); free(row); free(L.i); free(L.x); free(U.i); free(U.x); return false; }
    for (int64_t j = 0; j < n; j++) where[j] = -1;
    p->lp[0] = 0;
    p->up[0] = 0;
    bool ok = true;
    for (int64_t i = 0; i < n && ok; i++) {
        int64_t r0 = A->rp[i], r1 = A->rp[i + 1];
        bool has_diag = false;
        for (int64_t k = r0; k < r1; k++) {
            where[A->ci[k]] = k;
            vset(&p->s, row, A->ci[k], mget(A, k));
            if (A->ci[k] == i) has_diag = true;
        }
        if (!has_diag) { ok = false; break; }
        /* eliminate with the earlier rows, in column order (ci is sorted) */
        for (int64_t k = r0; k < r1; k++) {
            int64_t j = A->ci[k];
            if (j >= i) break;
            int64_t d = p->up[j];
            Cx piv = w == 1 ? cx(U.x[d], 0.0) : cx(U.x[2 * d], U.x[2 * d + 1]);
            Cx lij = cdiv(vget(&p->s, row, j), piv);
            vset(&p->s, row, j, lij);
            for (int64_t t = d + 1; t < p->up[j + 1]; t++) {
                int64_t c = U.i[t];
                if (where[c] < 0) continue;
                Cx u = w == 1 ? cx(U.x[t], 0.0) : cx(U.x[2 * t], U.x[2 * t + 1]);
                Cx cur = vget(&p->s, row, c);
                Cx pr = cmul(lij, u);
                vset(&p->s, row, c, cx(cur.re - pr.re, cur.im - pr.im));
            }
        }
        for (int64_t k = r0; k < r1 && ok; k++)
            if (A->ci[k] < i) ok = grow_push(&L, A->ci[k], vget(&p->s, row, A->ci[k]));
        Cx dg = vget(&p->s, row, i);
        if (!(cabs2(dg) > 0.0)) { ok = false; break; }
        ok = ok && grow_push(&U, i, dg);
        for (int64_t k = r0; k < r1 && ok; k++)
            if (A->ci[k] > i) ok = grow_push(&U, A->ci[k], vget(&p->s, row, A->ci[k]));
        p->lp[i + 1] = L.len;
        p->up[i + 1] = U.len;
        for (int64_t k = r0; k < r1; k++) {
            where[A->ci[k]] = -1;
            vset(&p->s, row, A->ci[k], cx(0.0, 0.0));
        }
    }
    free(where);
    free(row);
    p->li = L.i; p->lx = L.x;
    p->ui = U.i; p->ux = U.x;
    /* up[i] is the row start and the diagonal sits first in each U row:
     * prec_apply reads the diagonal at up[i]; the row ends at up[i + 1] */
    return ok;
}

/* ILUT(tau, fill) and ILUTP (column pivoting, permtol 0.5), after Saad,
 * "Iterative Methods for Sparse Linear Systems", sec. 10.4. */
static int cmp_mag_desc(const void* a, const void* b) {
    const double* x = a; const double* y = b;
    return (x[0] < y[0]) - (x[0] > y[0]);
}

static bool prec_ilut(const SpMat* A, Prec* p, double tau, int fill, bool pivot) {
    int64_t n = A->n;
    int w = p->s.w;
    Grow L = { 0, 0, 0, 0, 0, w }, U = { 0, 0, 0, 0, 0, w };
    p->lp = malloc(sizeof(int64_t) * (size_t)(n + 1));
    p->up = malloc(sizeof(int64_t) * (size_t)(n + 1));
    double* row = vnew(&p->s);
    int64_t* nzl = malloc(sizeof(int64_t) * (size_t)(n > 0 ? n : 1));   /* columns present */
    int64_t* mark = malloc(sizeof(int64_t) * (size_t)(n > 0 ? n : 1));
    double* sortbuf = malloc(sizeof(double) * 2 * (size_t)(n > 0 ? n : 1));
    int64_t* perm = NULL, *iperm = NULL;   /* perm[j]: A column at factor column j */
    if (pivot) {
        perm = malloc(sizeof(int64_t) * (size_t)(n > 0 ? n : 1));
        iperm = malloc(sizeof(int64_t) * (size_t)(n > 0 ? n : 1));
    }
    bool ok = p->lp && p->up && row && nzl && mark && sortbuf && (!pivot || (perm && iperm));
    if (ok) {
        for (int64_t j = 0; j < n; j++) mark[j] = -1;
        if (pivot) for (int64_t j = 0; j < n; j++) perm[j] = iperm[j] = j;
        p->lp[0] = 0;
        p->up[0] = 0;
    }
    for (int64_t i = 0; i < n && ok; i++) {
        /* scatter row i (in factor columns) */
        int64_t cnt = 0;
        double rnorm = 0.0;
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
            int64_t c = pivot ? iperm[A->ci[k]] : A->ci[k];
            Cx v = mget(A, k);
            vset(&p->s, row, c, v);
            mark[c] = i;
            nzl[cnt++] = c;
            rnorm += cabs2(v);
        }
        rnorm = sqrt(rnorm / (double)(A->rp[i + 1] - A->rp[i] > 0 ? A->rp[i + 1] - A->rp[i] : 1));
        double drop = tau * rnorm;
        /* eliminate the lower part in increasing column order; new fill can
         * appear below the current column, so pick the smallest each time */
        int64_t done = 0;
        int64_t* lower = malloc(sizeof(int64_t) * (size_t)(n > 0 ? n : 1));
        int64_t nlow = 0;
        if (!lower) { ok = false; break; }
        for (;;) {
            int64_t jmin = -1, at = -1;
            for (int64_t t = done; t < cnt; t++)
                if (nzl[t] < i && (jmin < 0 || nzl[t] < jmin)) { jmin = nzl[t]; at = t; }
            if (jmin < 0) break;
            /* move it to the done prefix */
            int64_t tmpc = nzl[done]; nzl[done] = nzl[at]; nzl[at] = tmpc;
            done++;
            Cx v = vget(&p->s, row, jmin);
            int64_t d = p->up[jmin];
            Cx piv = w == 1 ? cx(U.x[d], 0.0) : cx(U.x[2 * d], U.x[2 * d + 1]);
            Cx lij = cdiv(v, piv);
            if (sqrt(cabs2(lij)) < drop) { vset(&p->s, row, jmin, cx(0.0, 0.0)); continue; }
            vset(&p->s, row, jmin, lij);
            lower[nlow++] = jmin;
            for (int64_t t = d + 1; t < p->up[jmin + 1]; t++) {
                int64_t c = U.i[t];
                Cx u = w == 1 ? cx(U.x[t], 0.0) : cx(U.x[2 * t], U.x[2 * t + 1]);
                Cx pr = cmul(lij, u);
                if (mark[c] != i) {
                    mark[c] = i;
                    nzl[cnt++] = c;
                    vset(&p->s, row, c, cx(-pr.re, -pr.im));
                } else {
                    Cx cur = vget(&p->s, row, c);
                    vset(&p->s, row, c, cx(cur.re - pr.re, cur.im - pr.im));
                }
            }
        }
        /* keep the `fill` largest of L (beyond A's own count) */
        int64_t keepL = (A->rp[i + 1] - A->rp[i]) + fill;
        int64_t m = 0;
        for (int64_t t = 0; t < nlow; t++) {
            Cx v = vget(&p->s, row, lower[t]);
            sortbuf[2 * m] = sqrt(cabs2(v));
            sortbuf[2 * m + 1] = (double)lower[t];
            m++;
        }
        if (m > keepL) qsort(sortbuf, (size_t)m, 2 * sizeof(double), cmp_mag_desc);
        int64_t mk = m > keepL ? keepL : m;
        /* L row in increasing column order for the forward solve */
        for (int64_t a = 0; a < mk; a++)
            for (int64_t b = a + 1; b < mk; b++)
                if (sortbuf[2 * b + 1] < sortbuf[2 * a + 1]) {
                    double t0 = sortbuf[2 * a], t1 = sortbuf[2 * a + 1];
                    sortbuf[2 * a] = sortbuf[2 * b]; sortbuf[2 * a + 1] = sortbuf[2 * b + 1];
                    sortbuf[2 * b] = t0; sortbuf[2 * b + 1] = t1;
                }
        for (int64_t a = 0; a < mk && ok; a++) {
            int64_t c = (int64_t)sortbuf[2 * a + 1];
            ok = grow_push(&L, c, vget(&p->s, row, c));
        }
        free(lower);
        /* the upper part: diagonal (after optional pivoting) then the `fill`
         * largest off-diagonals */
        m = 0;
        for (int64_t t = 0; t < cnt; t++) {
            int64_t c = nzl[t];
            if (c < i) continue;
            Cx v = vget(&p->s, row, c);
            double mag = sqrt(cabs2(v));
            if (c != i && mag < drop) continue;
            sortbuf[2 * m] = mag;
            sortbuf[2 * m + 1] = (double)c;
            m++;
        }
        if (pivot) {
            /* swap in the largest entry when the diagonal is too small */
            double dmag = 0.0, best = 0.0;
            int64_t bc = -1;
            for (int64_t a = 0; a < m; a++) {
                int64_t c = (int64_t)sortbuf[2 * a + 1];
                if (c == i) dmag = sortbuf[2 * a];
                if (sortbuf[2 * a] > best) { best = sortbuf[2 * a]; bc = c; }
            }
            if (bc >= 0 && bc != i && 0.5 * best > dmag) {
                /* exchange factor columns i and bc everywhere from now on */
                Cx vi = vget(&p->s, row, i), vb = vget(&p->s, row, bc);
                vset(&p->s, row, i, vb);
                vset(&p->s, row, bc, vi);
                bool had_i = mark[i] == i;
                mark[i] = i;
                if (!had_i) { mark[bc] = i; }
                for (int64_t a = 0; a < m; a++) {
                    int64_t c = (int64_t)sortbuf[2 * a + 1];
                    if (c == i) sortbuf[2 * a + 1] = (double)bc;
                    else if (c == bc) sortbuf[2 * a + 1] = (double)i;
                }
                int64_t ai = perm[i], ab = perm[bc];
                perm[i] = ab; perm[bc] = ai;
                iperm[ab] = i; iperm[ai] = bc;
                /* earlier U rows refer to factor columns: swap i and bc there */
                for (int64_t t = 0; t < U.len; t++) {
                    if (U.i[t] == i) U.i[t] = bc;
                    else if (U.i[t] == bc) U.i[t] = i;
                }
            }
        }
        Cx dg = vget(&p->s, row, i);
        if (!(cabs2(dg) > 0.0)) dg = cx(drop > 0.0 ? drop : 1e-14, 0.0);   /* standard ILUT pivot fix */
        ok = ok && grow_push(&U, i, dg);
        int64_t keepU = (A->rp[i + 1] - A->rp[i]) + fill;
        int64_t mu = 0;
        for (int64_t a = 0; a < m; a++) {
            if ((int64_t)sortbuf[2 * a + 1] == i) continue;
            sortbuf[2 * mu] = sortbuf[2 * a];
            sortbuf[2 * mu + 1] = sortbuf[2 * a + 1];
            mu++;
        }
        if (mu > keepU) qsort(sortbuf, (size_t)mu, 2 * sizeof(double), cmp_mag_desc);
        for (int64_t a = 0; a < mu && a < keepU && ok; a++) {
            int64_t c = (int64_t)sortbuf[2 * a + 1];
            ok = grow_push(&U, c, vget(&p->s, row, c));
        }
        p->lp[i + 1] = L.len;
        p->up[i + 1] = U.len;
        for (int64_t t = 0; t < cnt; t++) vset(&p->s, row, nzl[t], cx(0.0, 0.0));
        vset(&p->s, row, i, cx(0.0, 0.0));
    }
    free(row); free(nzl); free(mark); free(sortbuf); free(iperm);
    p->li = L.i; p->lx = L.x;
    p->ui = U.i; p->ux = U.x;
    p->perm = perm;
    if (!ok) return false;
    /* the triangular solve visits U's off-diagonals as columns > i of the
     * factor; after pivoting every U entry column is a factor column */
    return true;
}

static bool prec_build(const SpMat* A, const KryOpts* o, KryPrec kind, Prec* p) {
    memset(p, 0, sizeof *p);
    p->s.n = A->n;
    p->s.w = A->cplx ? 2 : 1;
    p->kind = kind;
    p->tmp = vnew(&p->s);
    if (!p->tmp) return false;
    switch (kind) {
        case KPC_NONE: return true;
        case KPC_USER: p->user = o->user; p->user_ctx = o->user_ctx; return true;
        case KPC_JACOBI: return prec_jacobi(A, p);
        case KPC_ILU0: return prec_ilu0(A, p);
        case KPC_ILUT: return prec_ilut(A, p, o->ilut_droptol, o->ilut_fill, false);
        case KPC_ILUTP: return prec_ilut(A, p, o->ilut_droptol, o->ilut_fill, true);
        default: return false;
    }
}

/* ======================================================================== */
/* The iterations                                                           */
/* ======================================================================== */

static KryStatus run_cg(const SpMat* A, const double* b, double* x, const Prec* M,
                        double tol, int64_t maxit, int64_t* iters) {
    Vsp s = M->s;
    double *r = vnew(&s), *z = vnew(&s), *pv = vnew(&s), *q = vnew(&s);
    KryStatus st = KRY_NOMEM;
    if (!r || !z || !pv || !q) goto out;
    double bn = vnorm(&s, b);
    if (bn == 0.0) { memset(x, 0, sizeof(double) * (size_t)s.w * (size_t)s.n); *iters = 0; st = KRY_CONVERGED; goto out; }
    spla_spmv(A, x, r, false);
    for (int64_t i = 0; i < s.w * s.n; i++) r[i] = b[i] - r[i];
    if (!prec_apply(M, r, z)) { st = KRY_PRECFAIL; goto out; }
    vcopy(&s, z, pv);
    Cx rz = vdot(&s, r, z);
    st = KRY_MAXIT;
    for (int64_t it = 1; it <= maxit; it++) {
        *iters = it;
        spla_spmv(A, pv, q, false);
        Cx pq = vdot(&s, pv, q);
        if (cabs2(pq) == 0.0) { st = KRY_BREAKDOWN; break; }
        Cx alpha = cdiv(rz, pq);
        vaxpy(&s, alpha, pv, x);
        vaxpy(&s, cx(-alpha.re, -alpha.im), q, r);
        if (vnorm(&s, r) <= tol * bn) { st = KRY_CONVERGED; break; }
        if (!prec_apply(M, r, z)) { st = KRY_PRECFAIL; break; }
        Cx rz2 = vdot(&s, r, z);
        Cx beta = cdiv(rz2, rz);
        rz = rz2;
        /* p = z + beta p */
        vscale(&s, beta, pv);
        vaxpy(&s, cx(1.0, 0.0), z, pv);
    }
out:
    free(r); free(z); free(pv); free(q);
    return st;
}

static KryStatus run_bicgstab(const SpMat* A, const double* b, double* x, const Prec* M,
                              double tol, int64_t maxit, int64_t* iters) {
    Vsp s = M->s;
    double *r = vnew(&s), *r0 = vnew(&s), *pv = vnew(&s), *v = vnew(&s);
    double *ph = vnew(&s), *sv = vnew(&s), *sh = vnew(&s), *t = vnew(&s);
    KryStatus st = KRY_NOMEM;
    if (!r || !r0 || !pv || !v || !ph || !sv || !sh || !t) goto out;
    double bn = vnorm(&s, b);
    *iters = 0;
    if (bn == 0.0) { memset(x, 0, sizeof(double) * (size_t)s.w * (size_t)s.n); st = KRY_CONVERGED; goto out; }
    spla_spmv(A, x, r, false);
    for (int64_t i = 0; i < s.w * s.n; i++) r[i] = b[i] - r[i];
    if (vnorm(&s, r) <= tol * bn) { st = KRY_CONVERGED; goto out; }
    vcopy(&s, r, r0);
    Cx rho = cx(1.0, 0.0), alpha = cx(1.0, 0.0), omega = cx(1.0, 0.0);
    st = KRY_MAXIT;
    for (int64_t it = 1; it <= maxit; it++) {
        *iters = it;
        Cx rho1 = vdot(&s, r0, r);
        if (cabs2(rho1) == 0.0) { st = KRY_BREAKDOWN; break; }
        if (it == 1) {
            vcopy(&s, r, pv);
        } else {
            Cx beta = cmul(cdiv(rho1, rho), cdiv(alpha, omega));
            /* p = r + beta (p - omega v) */
            vaxpy(&s, cx(-omega.re, -omega.im), v, pv);
            vscale(&s, beta, pv);
            vaxpy(&s, cx(1.0, 0.0), r, pv);
        }
        rho = rho1;
        if (!prec_apply(M, pv, ph)) { st = KRY_PRECFAIL; break; }
        spla_spmv(A, ph, v, false);
        Cx r0v = vdot(&s, r0, v);
        if (cabs2(r0v) == 0.0) { st = KRY_BREAKDOWN; break; }
        alpha = cdiv(rho, r0v);
        vcopy(&s, r, sv);
        vaxpy(&s, cx(-alpha.re, -alpha.im), v, sv);
        if (vnorm(&s, sv) <= tol * bn) {
            vaxpy(&s, alpha, ph, x);
            st = KRY_CONVERGED;
            break;
        }
        if (!prec_apply(M, sv, sh)) { st = KRY_PRECFAIL; break; }
        spla_spmv(A, sh, t, false);
        Cx tt = vdot(&s, t, t);
        if (cabs2(tt) == 0.0) { st = KRY_BREAKDOWN; break; }
        omega = cdiv(vdot(&s, t, sv), tt);
        vaxpy(&s, alpha, ph, x);
        vaxpy(&s, omega, sh, x);
        vcopy(&s, sv, r);
        vaxpy(&s, cx(-omega.re, -omega.im), t, r);
        if (vnorm(&s, r) <= tol * bn) { st = KRY_CONVERGED; break; }
        if (cabs2(omega) == 0.0) { st = KRY_BREAKDOWN; break; }
    }
out:
    free(r); free(r0); free(pv); free(v); free(ph); free(sv); free(sh); free(t);
    return st;
}

static KryStatus run_gmres(const SpMat* A, const double* b, double* x, const Prec* M,
                           double tol, int64_t maxit, int restart, int64_t* iters) {
    Vsp s = M->s;
    int m = restart;
    size_t vlen = (size_t)s.w * (size_t)(s.n > 0 ? s.n : 1);
    double* V = calloc(vlen * (size_t)(m + 1), sizeof(double));     /* Krylov basis */
    double* Z = calloc(vlen * (size_t)m, sizeof(double));           /* preconditioned */
    Cx* H = calloc((size_t)(m + 1) * (size_t)m, sizeof(Cx));
    Cx* cs = calloc((size_t)m, sizeof(Cx));
    Cx* sn = calloc((size_t)m, sizeof(Cx));
    Cx* g = calloc((size_t)(m + 1), sizeof(Cx));
    Cx* yv = calloc((size_t)m, sizeof(Cx));
    double* w = vnew(&s);
    KryStatus st = KRY_NOMEM;
    *iters = 0;
    if (!V || !Z || !H || !cs || !sn || !g || !yv || !w) goto out;
    double bn = vnorm(&s, b);
    if (bn == 0.0) { memset(x, 0, sizeof(double) * vlen); st = KRY_CONVERGED; goto out; }
#define HH(i, j) H[(size_t)(i) * (size_t)m + (size_t)(j)]
    st = KRY_MAXIT;
    int64_t total = 0;
    while (total < maxit) {
        spla_spmv(A, x, w, false);
        for (int64_t i = 0; i < s.w * s.n; i++) w[i] = b[i] - w[i];
        double beta = vnorm(&s, w);
        if (beta <= tol * bn) { st = KRY_CONVERGED; break; }
        memcpy(V, w, sizeof(double) * vlen);
        vscale(&s, cx(1.0 / beta, 0.0), V);
        for (int i = 0; i <= m; i++) g[i] = cx(0.0, 0.0);
        g[0] = cx(beta, 0.0);
        int k = 0;
        bool conv = false;
        for (; k < m && total < maxit; k++) {
            total++;
            *iters = total;
            double* vk = V + vlen * (size_t)k;
            double* zk = Z + vlen * (size_t)k;
            if (!prec_apply(M, vk, zk)) { st = KRY_PRECFAIL; goto out; }
            spla_spmv(A, zk, w, false);
            /* modified Gram-Schmidt */
            for (int j = 0; j <= k; j++) {
                double* vj = V + vlen * (size_t)j;
                Cx h = vdot(&s, vj, w);
                HH(j, k) = h;
                vaxpy(&s, cx(-h.re, -h.im), vj, w);
            }
            double hn = vnorm(&s, w);
            HH(k + 1, k) = cx(hn, 0.0);
            if (hn > 0.0) {
                memcpy(V + vlen * (size_t)(k + 1), w, sizeof(double) * vlen);
                vscale(&s, cx(1.0 / hn, 0.0), V + vlen * (size_t)(k + 1));
            }
            /* apply the earlier rotations, then a new one zeroing H(k+1,k) */
            for (int j = 0; j < k; j++) {
                Cx a = HH(j, k), bb = HH(j + 1, k);
                Cx t1 = cmul(cconj(cs[j]), a), t2 = cmul(cconj(sn[j]), bb);
                Cx na = cx(t1.re + t2.re, t1.im + t2.im);
                Cx u1 = cmul(sn[j], a), u2 = cmul(cs[j], bb);
                HH(j + 1, k) = cx(-u1.re + u2.re, -u1.im + u2.im);
                HH(j, k) = na;
            }
            Cx a = HH(k, k), bb = HH(k + 1, k);
            double den = sqrt(cabs2(a) + cabs2(bb));
            if (den == 0.0) { cs[k] = cx(1.0, 0.0); sn[k] = cx(0.0, 0.0); }
            else { cs[k] = cx(a.re / den, a.im / den); sn[k] = cx(bb.re / den, bb.im / den); }
            {
                Cx t1 = cmul(cconj(cs[k]), a), t2 = cmul(cconj(sn[k]), bb);
                HH(k, k) = cx(t1.re + t2.re, t1.im + t2.im);
                HH(k + 1, k) = cx(0.0, 0.0);
                Cx gk = g[k];
                g[k] = cmul(cconj(cs[k]), gk);
                Cx u = cmul(sn[k], gk);
                g[k + 1] = cx(-u.re, -u.im);
            }
            if (sqrt(cabs2(g[k + 1])) <= tol * bn || hn == 0.0) { k++; conv = true; break; }
        }
        /* y = H(0:k,0:k)^-1 g, x += Z y */
        for (int i = k; i-- > 0;) {
            Cx acc = g[i];
            for (int j = i + 1; j < k; j++) {
                Cx pr = cmul(HH(i, j), yv[j]);
                acc = cx(acc.re - pr.re, acc.im - pr.im);
            }
            if (cabs2(HH(i, i)) == 0.0) { st = KRY_BREAKDOWN; goto out; }
            yv[i] = cdiv(acc, HH(i, i));
        }
        for (int j = 0; j < k; j++) vaxpy(&s, yv[j], Z + vlen * (size_t)j, x);
        if (conv) {
            /* confirm with the true residual (rounding in the recurrence) */
            spla_spmv(A, x, w, false);
            for (int64_t i = 0; i < s.w * s.n; i++) w[i] = b[i] - w[i];
            if (vnorm(&s, w) <= 10.0 * tol * bn) { st = KRY_CONVERGED; break; }
        }
    }
#undef HH
out:
    free(V); free(Z); free(H); free(cs); free(sn); free(g); free(yv); free(w);
    return st;
}

void kry_defaults(KryOpts* o) {
    memset(o, 0, sizeof *o);
    o->method = KRY_AUTO;
    o->prec = KPC_AUTO;
    o->tol = 1e-12;
    o->max_iter = 0;
    o->restart = 100;
    o->ilut_droptol = 1e-4;
    o->ilut_fill = 10;
}

KryStatus kry_solve(const SpMat* A, const double* b, double* x, const KryOpts* o, int64_t* iters) {
    *iters = 0;
    if (A->m != A->n) return KRY_BREAKDOWN;
    KryMethod meth = o->method == KRY_AUTO ? KRY_BICGSTAB : o->method;
    KryPrec kind = o->prec == KPC_AUTO ? KPC_ILU0 : o->prec;
    int64_t maxit = o->max_iter > 0 ? o->max_iter : (A->n < 100 ? 1000 : 10 * A->n);
    double tol = o->tol > 0.0 ? o->tol : 1e-12;
    Prec M;
    if (!prec_build(A, o, kind, &M)) {
        prec_free(&M);
        /* ILU0 needs every diagonal entry: fall back to ILUT, then none */
        if (kind == KPC_ILU0 && o->prec == KPC_AUTO) {
            if (!prec_build(A, o, KPC_ILUT, &M)) {
                prec_free(&M);
                if (!prec_build(A, o, KPC_NONE, &M)) { prec_free(&M); return KRY_NOMEM; }
            }
        } else {
            return KRY_PRECFAIL;
        }
    }
    KryStatus st;
    switch (meth) {
        case KRY_CG:    st = run_cg(A, b, x, &M, tol, maxit, iters); break;
        case KRY_GMRES: st = run_gmres(A, b, x, &M, tol, maxit, o->restart > 0 ? o->restart : 100, iters); break;
        default:        st = run_bicgstab(A, b, x, &M, tol, maxit, iters); break;
    }
    prec_free(&M);
    return st;
}

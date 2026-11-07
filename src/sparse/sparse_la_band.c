/* Mathilda -- sparse linear algebra: the band backend of SpFactor.
 *
 * A sparse direct solver for the builds without SuiteSparse (UMFPACK /
 * CHOLMOD), where the factor otherwise had only a dense LAPACK LU of the
 * densified matrix (O(n^3), and impossible beyond the densify cap -- the
 * Krylov fallback then answered).  Most sparse systems met in practice are
 * banded, or become banded under a symmetric renumbering: tridiagonal and
 * pentadiagonal matrices, finite-difference Laplacians on 2-D and 3-D grids
 * (bandwidth m and m^2 for an m-point side), block-tridiagonal systems.
 *
 *   1. Ordering.  The Reverse Cuthill-McKee renumbering of the symmetric
 *      pattern of A + A^T (George & Liu, Computer Solution of Large Sparse
 *      Positive Definite Systems, 1981, ch. 4): a breadth-first sweep from a
 *      pseudo-peripheral node of each connected component, neighbours taken
 *      in increasing degree, the order reversed.  The natural order is kept
 *      when it is no wider.
 *   2. Cost model.  The band factor costs ~n kl (kl + ku) multiply-adds in
 *      (2 kl + ku + 1) n doubles (LU) or ~n kd^2 / 2 in (kd + 1) n (Cholesky);
 *      the band is taken when it fits the densify cap and costs at most an
 *      eighth of the dense LU (n^3 / 3), else the caller's dense path runs.
 *   3. Factor.  LAPACK dpbtrf (zpbtrf) for a Hermitian matrix with a
 *      positive diagonal, falling back to dgbtrf (zgbtrf), the band LU with
 *      partial pivoting, when it is not positive definite.  Both are
 *      backward stable: the answer is as accurate as the dense LU's.
 *
 * The symmetric permutation P leaves the determinant and every reciprocal
 * condition estimate unchanged: A x = b is (P A P^T)(P x) = P b.
 */
#include "sparse_la.h"
#include "sparse_la_band.h"
#include "lapack.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* The densify cap (elements), the same as the dense backend's. */
#define SPB_STORE_MAX ((int64_t)1 << 27)

struct SpBand {
    int64_t n;
    bool cplx;
    bool chol;          /* pbtrf factor (lower band, ldab = kd + 1)         */
    int kl, ku, ldab;   /* gbtrf: ldab = 2 kl + ku + 1; pbtrf: kl = ku = kd */
    double* ab;
    int* ipiv;
    int64_t* perm;      /* new index -> old index (NULL: natural order)     */
    double rcond;
};

void spb_free(SpBand* b) {
    if (!b) return;
    free(b->ab);
    free(b->ipiv);
    free(b->perm);
    free(b);
}

bool spb_is_chol(const SpBand* b) { return b->chol; }
double spb_rcond(const SpBand* b) { return b->rcond; }

/* ---- the symmetric pattern and RCM ---------------------------------------- */

typedef struct {
    int64_t n;
    int64_t* p;         /* n + 1 */
    int64_t* j;         /* p[n] neighbours, the diagonal excluded */
} Graph;

static void graph_free(Graph* g) { free(g->p); free(g->j); }

/* The adjacency of A + A^T (duplicates removed). */
static bool graph_build(const SpMat* A, Graph* g) {
    int64_t n = A->n;
    g->n = n;
    g->p = calloc((size_t)n + 1, sizeof(int64_t));
    int64_t* cnt = calloc((size_t)n + 1, sizeof(int64_t));
    if (!g->p || !cnt) { free(cnt); free(g->p); g->p = NULL; g->j = NULL; return false; }
    for (int64_t i = 0; i < n; i++)
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
            int64_t c = A->ci[k];
            if (c == i) continue;
            cnt[i]++;
            cnt[c]++;
        }
    int64_t tot = 0;
    for (int64_t i = 0; i < n; i++) { g->p[i] = tot; tot += cnt[i]; }
    g->p[n] = tot;
    g->j = malloc(sizeof(int64_t) * (size_t)(tot > 0 ? tot : 1));
    if (!g->j) { free(cnt); free(g->p); g->p = NULL; return false; }
    for (int64_t i = 0; i < n; i++) cnt[i] = g->p[i];
    for (int64_t i = 0; i < n; i++)
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
            int64_t c = A->ci[k];
            if (c == i) continue;
            g->j[cnt[i]++] = c;
            g->j[cnt[c]++] = i;
        }
    /* deduplicate each list in place with a marker array */
    int64_t* mark = cnt;              /* reuse: mark[v] = last row that saw v */
    for (int64_t i = 0; i < n; i++) mark[i] = -1;
    int64_t w = 0;
    for (int64_t i = 0; i < n; i++) {
        int64_t lo = g->p[i], hi = g->p[i + 1];
        g->p[i] = w;
        for (int64_t k = lo; k < hi; k++) {
            int64_t v = g->j[k];
            if (mark[v] == i) continue;
            mark[v] = i;
            g->j[w++] = v;
        }
    }
    g->p[n] = w;
    free(cnt);
    return true;
}

static int64_t deg(const Graph* g, int64_t v) { return g->p[v + 1] - g->p[v]; }

/* Breadth-first levels from `root`: fills `order` (the visit order, `*cnt`
 * vertices, `*last_lo` the first of the last level) and returns the
 * eccentricity (the number of levels - 1).  `lvl` is scratch, -1 = unseen;
 * the caller resets the entries of `order` to -1 afterwards. */
static int64_t bfs_levels(const Graph* g, int64_t root, int64_t* order, int64_t* lvl,
                          int64_t* cnt, int64_t* last_lo) {
    int64_t head = 0, tail = 0;
    order[tail++] = root;
    lvl[root] = 0;
    int64_t ecc = 0;
    *last_lo = 0;
    while (head < tail) {
        int64_t v = order[head++];
        if (lvl[v] > ecc) { ecc = lvl[v]; *last_lo = head - 1; }
        for (int64_t k = g->p[v]; k < g->p[v + 1]; k++) {
            int64_t u = g->j[k];
            if (lvl[u] != -1) continue;
            lvl[u] = lvl[v] + 1;
            order[tail++] = u;
        }
    }
    *cnt = tail;
    return ecc;
}

/* qsort context: the graph whose degrees order the neighbours (the band
 * factor runs on the main thread only, never re-entrantly) */
static const Graph* g_sort_graph;
static int cmp_deg_ctx_cmp(const void* a, const void* b) {
    int64_t x = *(const int64_t*)a, y = *(const int64_t*)b;
    int64_t dx = deg(g_sort_graph, x), dy = deg(g_sort_graph, y);
    if (dx != dy) return dx < dy ? -1 : 1;
    return x < y ? -1 : (x > y);
}

/* RCM: perm[new] = old.  False on allocation failure. */
static bool rcm_order(const Graph* g, int64_t* perm) {
    int64_t n = g->n;
    int64_t* lvl = malloc(sizeof(int64_t) * (size_t)(n > 0 ? n : 1));
    int64_t* tmp = malloc(sizeof(int64_t) * (size_t)(n > 0 ? n : 1));
    char* done = calloc((size_t)(n > 0 ? n : 1), 1);
    if (!lvl || !tmp || !done) { free(lvl); free(tmp); free(done); return false; }
    for (int64_t i = 0; i < n; i++) lvl[i] = -1;
    int64_t pos = 0;
    g_sort_graph = g;
    for (int64_t s = 0; s < n; s++) {
        if (done[s]) continue;
        /* a pseudo-peripheral root of s's component (George-Liu): restart the
         * sweep from the lowest-degree vertex of the last level while the
         * eccentricity grows.  Each sweep resets the levels it set, so the
         * whole ordering is O(nnz) per sweep, not O(n) per component. */
        int64_t root = s, cnt = 0, last_lo = 0;
        int64_t ecc = bfs_levels(g, root, tmp, lvl, &cnt, &last_lo);
        for (int64_t k = 0; k < cnt; k++) lvl[tmp[k]] = -1;
        for (int it = 0; it < 8 && ecc > 0; it++) {
            int64_t best = tmp[last_lo];
            for (int64_t k = last_lo; k < cnt; k++)
                if (deg(g, tmp[k]) < deg(g, best)) best = tmp[k];
            int64_t cnt2 = 0, lo2 = 0;
            int64_t e2 = bfs_levels(g, best, tmp, lvl, &cnt2, &lo2);
            for (int64_t k = 0; k < cnt2; k++) lvl[tmp[k]] = -1;
            cnt = cnt2;
            last_lo = lo2;
            if (e2 <= ecc) break;
            root = best;
            ecc = e2;
        }
        /* Cuthill-McKee from root: neighbours in increasing degree */
        int64_t head = pos, tail = pos;
        perm[tail++] = root;
        done[root] = 1;
        while (head < tail) {
            int64_t v = perm[head++];
            int64_t start = tail;
            for (int64_t k = g->p[v]; k < g->p[v + 1]; k++) {
                int64_t u = g->j[k];
                if (done[u]) continue;
                done[u] = 1;
                perm[tail++] = u;
            }
            if (tail - start > 1)
                qsort(perm + start, (size_t)(tail - start), sizeof(int64_t), cmp_deg_ctx_cmp);
        }
        pos = tail;
    }
    g_sort_graph = NULL;
    /* reverse */
    for (int64_t i = 0, j = n - 1; i < j; i++, j--) { int64_t t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
    free(lvl);
    free(tmp);
    free(done);
    return pos == n;
}

/* Lower / upper bandwidth of A under the numbering q (old -> new; NULL =
 * identity). */
static void bandwidth(const SpMat* A, const int64_t* q, int64_t* kl, int64_t* ku) {
    int64_t lo = 0, up = 0;
    for (int64_t i = 0; i < A->n; i++) {
        int64_t ni = q ? q[i] : i;
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
            int64_t nj = q ? q[A->ci[k]] : A->ci[k];
            int64_t d = ni - nj;
            if (d > lo) lo = d;
            if (-d > up) up = -d;
        }
    }
    *kl = lo;
    *ku = up;
}

/* ---- factor ----------------------------------------------------------------- */

/* Fill the band storage: LU layout (row kl + ku + i - j) or the lower
 * Cholesky layout (row i - j, i >= j). */
static void band_fill(const SpMat* A, const int64_t* q, SpBand* b) {
    int w = b->cplx ? 2 : 1;
    for (int64_t i = 0; i < A->n; i++) {
        int64_t ni = q ? q[i] : i;
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
            int64_t nj = q ? q[A->ci[k]] : A->ci[k];
            int64_t row;
            if (b->chol) {
                if (ni < nj) continue;
                row = ni - nj;
            } else {
                row = (int64_t)b->kl + b->ku + ni - nj;
            }
            size_t at = (size_t)w * ((size_t)row + (size_t)nj * (size_t)b->ldab);
            b->ab[at] += A->x[w * k];
            if (b->cplx) b->ab[at + 1] += A->x[w * k + 1];
        }
    }
}

/* ---- narrow real bands, in-house ------------------------------------------
 *
 * For a band of a few diagonals LAPACK's unblocked kernels (?pbtf2, ?gbtf2,
 * ?tbsv) make several BLAS calls of length kl per column, and the call
 * overhead -- not the arithmetic -- is the cost: a 10000 x 10000
 * pentadiagonal Cholesky spent ~1 ms in 30000 calls.  Up to NARROW total
 * diagonals the same algorithms run as plain loops, in the same storage and
 * with the same pivot convention, so the factor is LAPACK's to the bit
 * pattern of operations (column-by-column, the same order of updates). */
#define NARROW 16

#define AB(r, c) ab[(size_t)(r) + (size_t)(c) * (size_t)ldab]

/* ?pbtf2, lower: A = L L^T.  Returns 0, or j + 1 when the leading minor of
 * order j + 1 is not positive. */
static int narrow_pbtrf(int n, int kd, double* ab, int ldab) {
    for (int j = 0; j < n; j++) {
        double ajj = AB(0, j);
        if (!(ajj > 0.0)) return j + 1;
        ajj = sqrt(ajj);
        AB(0, j) = ajj;
        int kn = kd < n - 1 - j ? kd : n - 1 - j;
        double r = 1.0 / ajj;
        for (int i = 1; i <= kn; i++) AB(i, j) *= r;
        /* trailing update: A(j+1+c.., j+1+c) -= l(c) l(r) */
        for (int c = 0; c < kn; c++) {
            double lc = AB(1 + c, j);
            for (int rr = c; rr < kn; rr++) AB(rr - c, j + 1 + c) -= AB(1 + rr, j) * lc;
        }
    }
    return 0;
}

static void narrow_pbtrs(int n, int kd, const double* ab, int ldab, double* x) {
    for (int j = 0; j < n; j++) {               /* L y = b */
        double v = x[j] / AB(0, j);
        x[j] = v;
        int kn = kd < n - 1 - j ? kd : n - 1 - j;
        for (int i = 1; i <= kn; i++) x[j + i] -= AB(i, j) * v;
    }
    for (int j = n - 1; j >= 0; j--) {          /* L^T x = y */
        int kn = kd < n - 1 - j ? kd : n - 1 - j;
        double v = x[j];
        for (int i = 1; i <= kn; i++) v -= AB(i, j) * x[j + i];
        x[j] = v / AB(0, j);
    }
}

/* ?gbtf2: P A = L U with partial pivoting, kl sub- and ku superdiagonals,
 * rows 0..kl-1 of the storage holding the fill.  ipiv 1-based as LAPACK's.
 * Returns 0, or j + 1 for the first exactly zero pivot. */
static int narrow_gbtrf(int n, int kl, int ku, double* ab, int ldab, int* ipiv) {
    int kv = ku + kl, info = 0, ju = 0;
    for (int j = 0; j < n; j++) {
        int km = kl < n - 1 - j ? kl : n - 1 - j;
        int jp = 0;
        double best = fabs(AB(kv, j));
        for (int i = 1; i <= km; i++)
            if (fabs(AB(kv + i, j)) > best) { best = fabs(AB(kv + i, j)); jp = i; }
        ipiv[j] = jp + j + 1;
        if (AB(kv + jp, j) != 0.0) {
            int last = j + ku + jp < n - 1 ? j + ku + jp : n - 1;
            if (last > ju) ju = last;
            if (jp != 0)                          /* swap rows j and j + jp, cols j..ju */
                for (int c = j; c <= ju; c++) {
                    double t = AB(kv + jp + j - c, c);
                    AB(kv + jp + j - c, c) = AB(kv + j - c, c);
                    AB(kv + j - c, c) = t;
                }
            if (km > 0) {
                double r = 1.0 / AB(kv, j);
                for (int i = 1; i <= km; i++) AB(kv + i, j) *= r;
                for (int c = j + 1; c <= ju; c++) {
                    double t = AB(kv + j - c, c);
                    if (t != 0.0)
                        for (int i = 1; i <= km; i++) AB(kv + i + j - c, c) -= AB(kv + i, j) * t;
                }
            }
        } else if (info == 0) {
            info = j + 1;
        }
    }
    return info;
}

/* ?gbtrs on one right-hand side: A x = b (trans false) or A^T x = b. */
static void narrow_gbtrs(bool trans, int n, int kl, int ku, const double* ab, int ldab,
                         const int* ipiv, double* x) {
    int kv = ku + kl;
    if (!trans) {
        for (int j = 0; j < n; j++) {             /* L y = P b */
            int km = kl < n - 1 - j ? kl : n - 1 - j;
            int p = ipiv[j] - 1;
            if (p != j) { double t = x[j]; x[j] = x[p]; x[p] = t; }
            double v = x[j];
            if (v != 0.0) for (int i = 1; i <= km; i++) x[j + i] -= AB(kv + i, j) * v;
        }
        for (int j = n - 1; j >= 0; j--) {        /* U x = y: U(j, c) at row kv + j - c */
            double v = x[j] / AB(kv, j);
            x[j] = v;
            int i0 = j - kv > 0 ? j - kv : 0;
            if (v != 0.0) for (int i = i0; i < j; i++) x[i] -= AB(kv + i - j, j) * v;
        }
    } else {
        for (int j = 0; j < n; j++) {             /* U^T y = b */
            int i0 = j - kv > 0 ? j - kv : 0;
            double v = x[j];
            for (int i = i0; i < j; i++) v -= AB(kv + i - j, j) * x[i];
            x[j] = v / AB(kv, j);
        }
        for (int j = n - 1; j >= 0; j--) {        /* L^T z = y, then P^T */
            int km = kl < n - 1 - j ? kl : n - 1 - j;
            double v = x[j];
            for (int i = 1; i <= km; i++) v -= AB(kv + i, j) * x[j + i];
            x[j] = v;
            int p = ipiv[j] - 1;
            if (p != j) { double t = x[j]; x[j] = x[p]; x[p] = t; }
        }
    }
}

#undef AB

static bool is_narrow(const SpBand* b) {
    return !b->cplx && (b->chol ? 2 * b->kl : b->kl + b->ku) <= NARROW;
}

static SpfStatus band_factor(const SpMat* A, const int64_t* q, SpBand* b) {
#ifndef USE_LAPACK
    /* No LAPACK in this build (the native Windows package): the in-house
     * narrow-band loops still factor a narrow real band; a wider or complex
     * one is SPF_TOOBIG, which every caller reads as "use the next method". */
    if (!is_narrow(b)) return SPF_TOOBIG;
#endif
    int w = b->cplx ? 2 : 1;
    size_t len = (size_t)w * (size_t)b->ldab * (size_t)b->n;
    free(b->ab);
    b->ab = calloc(len > 0 ? len : 1, sizeof(double));
    if (!b->ab) return SPF_NOMEM;
    band_fill(A, q, b);
    int n = (int)b->n;
    int info;
    bool narrow = is_narrow(b);
    if (b->chol) {
        info = narrow ? narrow_pbtrf(n, b->kl, b->ab, b->ldab)
             : b->cplx ? mat_lapack_zpbtrf(n, b->kl, b->ab, b->ldab)
                       : mat_lapack_dpbtrf(n, b->kl, b->ab, b->ldab);
        if (info < 0) return SPF_FAILED;
        if (info > 0) return SPF_NOTPD;
    } else {
        free(b->ipiv);
        b->ipiv = malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
        if (!b->ipiv) return SPF_NOMEM;
        info = narrow ? narrow_gbtrf(n, b->kl, b->ku, b->ab, b->ldab, b->ipiv)
             : b->cplx ? mat_lapack_zgbtrf(n, b->kl, b->ku, b->ab, b->ldab, b->ipiv)
                       : mat_lapack_dgbtrf(n, b->kl, b->ku, b->ab, b->ldab, b->ipiv);
        if (info < 0) return SPF_FAILED;
        if (info > 0) return SPF_SINGULAR;
    }
    /* the pivot-magnitude ratio min |u_jj| / max |u_jj| (squared for the
     * Cholesky factor), the estimate a sparse LU gives (UMFPACK's RCOND) --
     * and what the reference's sparse ::luc reads; a 1-norm estimator
     * (?gbcon) costs several band solves of BLAS-1 calls and would double the
     * time of a narrow band */
    size_t row = b->chol ? 0 : (size_t)(b->kl + b->ku);
    double lo = HUGE_VAL, hi = 0.0;
    for (int64_t j = 0; j < b->n; j++) {
        size_t at = (size_t)w * (row + (size_t)j * (size_t)b->ldab);
        double m = b->cplx ? hypot(b->ab[at], b->ab[at + 1]) : fabs(b->ab[at]);
        if (m < lo) lo = m;
        if (m > hi) hi = m;
    }
    double rc = hi > 0.0 ? lo / hi : 0.0;
    b->rcond = b->chol ? rc * rc : rc;
    if (!(b->rcond > 0.0)) return SPF_SINGULAR;
    return SPF_OK;
}

SpBand* spb_factor(const SpMat* A, SpbMode mode, SpfStatus* st) {
    *st = SPF_TOOBIG;
    int64_t n = A->n;
    if (A->m != n || n < 2 || n > 0x7fffffff) return NULL;

    /* the natural order, then RCM when that is wide */
    int64_t kl0, ku0;
    bandwidth(A, NULL, &kl0, &ku0);
    int64_t* perm = NULL;
    int64_t kl = kl0, ku = ku0;
    /* a narrow natural band already factors in O(n) -- the renumbering would
     * cost more than it could save */
    if (kl0 + ku0 > 16) {
        Graph g;
        if (graph_build(A, &g)) {
            int64_t* p = malloc(sizeof(int64_t) * (size_t)n);
            int64_t* q = malloc(sizeof(int64_t) * (size_t)n);
            if (p && q && rcm_order(&g, p)) {
                for (int64_t i = 0; i < n; i++) q[p[i]] = i;
                int64_t kl1, ku1;
                bandwidth(A, q, &kl1, &ku1);
                if (kl1 + ku1 < kl0 + ku0) { kl = kl1; ku = ku1; perm = p; p = NULL; }
            }
            free(p);
            free(q);
            graph_free(&g);
        }
    }

    bool want_chol = (mode == SPB_CHOL || mode == SPB_CHOL_LU) && kl == ku;
    double dn = (double)n;
    double dense_cost = dn * dn * dn / 3.0;
    double lu_cost = dn * (double)kl * (double)(kl + ku + 1);
    double lu_store = (double)(2 * kl + ku + 1) * dn;
    double ch_cost = dn * (double)kl * (double)(kl + 1) / 2.0;
    double ch_store = (double)(kl + 1) * dn;
    bool lu_ok = lu_store <= (double)SPB_STORE_MAX && lu_cost <= dense_cost / 8.0;
    bool ch_ok = want_chol && ch_store <= (double)SPB_STORE_MAX && ch_cost <= dense_cost / 8.0;
    if (mode == SPB_CHOL ? !ch_ok : !(lu_ok || ch_ok)) { free(perm); return NULL; }

    SpBand* b = calloc(1, sizeof *b);
    if (!b) { free(perm); *st = SPF_NOMEM; return NULL; }
    b->n = n;
    b->cplx = A->cplx;
    b->rcond = -1.0;
    int64_t* q = NULL;
    if (perm) {
        q = malloc(sizeof(int64_t) * (size_t)n);
        if (!q) { free(perm); free(b); *st = SPF_NOMEM; return NULL; }
        for (int64_t i = 0; i < n; i++) q[perm[i]] = i;
    }
    b->perm = perm;
    SpfStatus s = SPF_FAILED;
    if (ch_ok) {
        b->chol = true;
        b->kl = b->ku = (int)kl;
        b->ldab = (int)kl + 1;
        s = band_factor(A, q, b);
        if (s == SPF_OK || mode == SPB_CHOL || !lu_ok) {
            free(q);
            if (s != SPF_OK) { *st = s == SPF_NOTPD && mode != SPB_CHOL ? SPF_TOOBIG : s; spb_free(b); return NULL; }
            *st = SPF_OK;
            return b;
        }
    }
    b->chol = false;
    b->kl = (int)kl;
    b->ku = (int)ku;
    b->ldab = (int)(2 * kl + ku + 1);
    s = band_factor(A, q, b);
    free(q);
    if (s != SPF_OK) { *st = s; spb_free(b); return NULL; }
    *st = SPF_OK;
    return b;
}

/* ---- solve / det ------------------------------------------------------------- */

bool spb_solve(const SpBand* b, const double* B, double* X, int64_t nrhs, bool trans) {
    int w = b->cplx ? 2 : 1;
    int64_t n = b->n;
    size_t len = (size_t)w * (size_t)n;
    double* y = malloc(sizeof(double) * len * (size_t)(nrhs > 0 ? nrhs : 1));
    if (!y) return false;
    /* y = P b (and conj for the transposed Hermitian Cholesky system:
     * A^T = conj(A) there, so A^T x = b is A conj(x) = conj(b)) */
    bool conj = b->chol && b->cplx && trans;
    for (int64_t r = 0; r < nrhs; r++) {
        const double* src = B + (size_t)r * len;
        double* dst = y + (size_t)r * len;
        for (int64_t i = 0; i < n; i++) {
            int64_t o = b->perm ? b->perm[i] : i;
            dst[w * i] = src[w * o];
            if (w == 2) dst[w * i + 1] = conj ? -src[w * o + 1] : src[w * o + 1];
        }
    }
    int info = 0;
    if (is_narrow(b)) {
        for (int64_t r = 0; r < nrhs; r++) {
            double* yr = y + (size_t)r * len;
            if (b->chol) narrow_pbtrs((int)n, b->kl, b->ab, b->ldab, yr);
            else narrow_gbtrs(trans, (int)n, b->kl, b->ku, b->ab, b->ldab, b->ipiv, yr);
        }
    } else if (b->chol)
        info = b->cplx ? mat_lapack_zpbtrs((int)n, b->kl, (int)nrhs, b->ab, b->ldab, y, (int)n)
                       : mat_lapack_dpbtrs((int)n, b->kl, (int)nrhs, b->ab, b->ldab, y, (int)n);
    else
        info = b->cplx ? mat_lapack_zgbtrs(trans ? 'T' : 'N', (int)n, b->kl, b->ku, (int)nrhs, b->ab, b->ldab, b->ipiv, y, (int)n)
                       : mat_lapack_dgbtrs(trans ? 'T' : 'N', (int)n, b->kl, b->ku, (int)nrhs, b->ab, b->ldab, b->ipiv, y, (int)n);
    if (info != 0) { free(y); return false; }
    for (int64_t r = 0; r < nrhs; r++) {
        const double* src = y + (size_t)r * len;
        double* dst = X + (size_t)r * len;
        for (int64_t i = 0; i < n; i++) {
            int64_t o = b->perm ? b->perm[i] : i;
            dst[w * o] = src[w * i];
            if (w == 2) dst[w * o + 1] = conj ? -src[w * i + 1] : src[w * i + 1];
        }
    }
    free(y);
    return true;
}

bool spb_det(const SpBand* b, double* re, double* im) {
    int w = b->cplx ? 2 : 1;
    double dr = 1.0, di = 0.0;
    for (int64_t j = 0; j < b->n; j++) {
        size_t row = b->chol ? 0 : (size_t)(b->kl + b->ku);
        size_t at = (size_t)w * (row + (size_t)j * (size_t)b->ldab);
        double ar = b->ab[at], ai = b->cplx ? b->ab[at + 1] : 0.0;
        if (b->chol) {
            /* det = prod l_jj^2, the diagonal real and positive */
            dr *= ar * ar;
            continue;
        }
        if (b->ipiv[j] - 1 != j) { dr = -dr; di = -di; }
        double t = dr * ar - di * ai;
        di = dr * ai + di * ar;
        dr = t;
    }
    *re = dr;
    *im = b->cplx ? di : 0.0;
    return isfinite(dr) && isfinite(di);
}

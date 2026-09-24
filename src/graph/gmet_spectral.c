/* gmet_spectral.c - spectral / random-walk centralities.
 *
 *   PageRankCentrality[g] / [g, a]
 *   EigenvectorCentrality[g] / [g, "In"] / [g, "Out"]
 *   KatzCentrality[g, a] / [g, a, b]
 *   HITSCentrality[g]
 *
 * Every one ignores EdgeWeight (Wolfram gives identical answers with and
 * without it) and treats an undirected edge as an arc in each direction.
 * Results are machine reals (packed), matching Mathematica 15, whose values
 * were used to pin down each convention:
 *
 *   PageRank: x = a P^T x + (1 - a)/n, where P is the row-stochastic walk
 *     matrix and a vertex with no outgoing arc (dangling) jumps to every
 *     vertex uniformly; sum(x) = 1. Default a = 0.85; 0 <= a <= 1. Power
 *     iteration to an L1 change below 1e-14 (Wolfram stops near 1e-9, so the
 *     two agree to about 9 digits).
 *   Eigenvector: computed per strongly connected component. Each non-trivial
 *     component C gets its Perron vector (left, x_v ~ sum over arcs u->v, for
 *     "In", the default; right for "Out"), normalized to sum (|C| - 1) /
 *     sum over components (|C'| - 1); vertices in trivial components (a single
 *     vertex) get 0. That weighting is not documented by Wolfram but
 *     reproduces its output exactly on every disconnected and non-strongly-
 *     connected case tried (e.g. K4+K3+K2 gives each clique total 3:2:1).
 *     Perron vectors come from restarted Arnoldi (see "Block Perron vectors"
 *     below); a block whose iteration does not converge leaves the call
 *     unevaluated rather than answered inaccurately.
 *   Katz: x = a A^T x + b, i.e. x_v = b_v + a sum_{u->v} x_u; b defaults to 1
 *     and may be a number or a list. Jacobi iteration when it converges
 *     (a * spectral radius < 1), else a dense LAPACK solve (n <= 4000), which
 *     like Wolfram answers even for a beyond the convergence radius (K4 with
 *     a = 1/2 gives -2). A singular system is left unevaluated (Wolfram:
 *     KatzCentrality::nosol). With no edges, or an exact a = 0, the answer is
 *     b itself, exactly, as in Wolfram.
 *   HITS: {h, a} with a = A h (not renormalized) and h the same block
 *     construction as Eigenvector applied to A^T A: its blocks are the classes
 *     of "share an in-neighbour" (connected components of the co-citation
 *     graph), each gets its Perron vector of A^T A scaled to total
 *     (|C| - 1) / sum (|C'| - 1), singletons get 0. Again undocumented, but it
 *     reproduces Wolfram exactly, including every degenerate case tried
 *     (PathGraph[5], StarGraph, disconnected and directed graphs) and the all-
 *     zero answer Wolfram gives for a mixed graph such as {1->2, 2<->3}, whose
 *     classes are all singletons.
 */

#include "graph_metrics.h"
#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "linalg/lapack.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* A real numeric argument as a double; 0 if not numeric. *exact_zero set when
 * the argument is the exact integer 0. */
static int read_real(const Expr* e, double* x, int* exact_zero) {
    if (exact_zero) *exact_zero = (e->type == EXPR_INTEGER && e->data.integer == 0);
    double d = graph_weight_to_double(e);
    if (isnan(d) || isinf(d)) return 0;
    *x = d;
    return 1;
}

/* ---- PageRank -------------------------------------------------------------- */

typedef struct {
    const GmetCSR* in;
    const double* y;      /* x[u] / outdeg[u] */
    double* xn;
    double alpha, base;
} PrCtx;

static void pr_pull(void* vctx, int tid, int64_t lo, int64_t hi) {
    (void)tid;
    PrCtx* c = (PrCtx*)vctx;
    const int64_t* off = c->in->off;
    const int* adj = c->in->adj;
    for (int64_t v = lo; v < hi; v++) {
        double s = 0.0;
        for (int64_t p = off[v]; p < off[v + 1]; p++) s += c->y[adj[p]];
        c->xn[v] = c->alpha * s + c->base;
    }
}

Expr* builtin_pagerank_centrality(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    double alpha = 0.85;
    if (argc == 2 && !read_real(res->data.function.args[1], &alpha, NULL)) return NULL;
    if (alpha < 0.0 || alpha > 1.0) return NULL;
    if (!graph_is_valid(g)) return NULL;
    Expr* cached = gmet_cache_get("PageRankCentrality", res);
    if (cached) return cached;

    GmetCSR in;
    if (!gmet_csr_build(g, GMET_IN, NULL, 0, &in)) return NULL;
    int n = in.n;
    size_t nn = (size_t)(n > 0 ? n : 1);
    double* outdeg = calloc(nn, sizeof(double));
    double* x = malloc(nn * sizeof(double));
    double* xn = malloc(nn * sizeof(double));
    double* y = malloc(nn * sizeof(double));
    if (!outdeg || !x || !xn || !y) {
        free(outdeg); free(x); free(xn); free(y); gmet_csr_free(&in); return NULL;
    }
    for (int64_t p = 0; p < in.narcs; p++) outdeg[in.adj[p]] += 1.0;
    for (int i = 0; i < n; i++) { x[i] = 1.0 / n; outdeg[i] = outdeg[i] > 0 ? 1.0 / outdeg[i] : 0.0; }
    PrCtx c = { &in, y, xn, alpha, 0.0 };
    int nt = gmet_thread_count((double)(in.narcs + n) * 20.0);
    for (int it = 0; it < 10000 && n > 0; it++) {
        double dang = 0.0;
        for (int i = 0; i < n; i++) {
            if (outdeg[i] == 0.0) dang += x[i];
            y[i] = x[i] * outdeg[i];
        }
        c.base = (alpha * dang + (1.0 - alpha)) / n;
        if (nt > 1) gmet_parallel_for(n, 4096, nt, pr_pull, &c);
        else pr_pull(&c, 0, 0, n);
        /* Renormalize BEFORE measuring the step: the teleport sum carries a
         * ~1e-13 rounding deficit on large graphs, which would otherwise
         * masquerade as a step that never shrinks. */
        double sum = 0.0;
        for (int i = 0; i < n; i++) sum += xn[i];
        double diff = 0.0;
        for (int i = 0; i < n; i++) { xn[i] /= sum; diff += fabs(xn[i] - x[i]); }
        double* t = x; x = xn; xn = t; c.xn = xn;
        if (diff < 1e-14) break;
    }
    Expr* out = gmet_real_vector(x, n);
    free(outdeg); free(x); free(xn); free(y);
    gmet_csr_free(&in);
    gmet_cache_put("PageRankCentrality", res, out);
    return out;
}

/* ---- Block Perron vectors ---------------------------------------------------
 * Both EigenvectorCentrality and HITSCentrality reduce to: split a non-negative
 * matrix into independent diagonal blocks, take each block's Perron vector,
 * and combine the blocks with weight (block size - 1). A block is given as a
 * local CSR over its k vertices; the operator is y = B x (one stage, row i of
 * B lists the columns summed into y_i) or y = B^T B x (two stage, for HITS'
 * A^T A, with B's rows the in-neighbours' out-lists).
 *
 * The Perron vector comes from explicitly restarted Arnoldi (Krylov dimension
 * up to 40, full re-orthogonalization by BLAS dgemv, LAPACK dgeev on the small Hessenberg
 * matrix, restart from the Ritz vector) until the residual ||M x - t x|| is
 * below 1e-13 |t|. Unlike plain power iteration this does not stall when the
 * second eigenvalue is close to the first (nearly bipartite graphs, long
 * paths), and the residual test, unlike a small-step test, bounds the error.
 * Blocks of at most 64 vertices are solved densely. Without LAPACK a shifted
 * power iteration with the same residual test is used, and a block that does
 * not converge makes the head return unevaluated. */

typedef struct {
    int k;
    int two_stage;
    int nrows;              /* rows of B (== k when one stage)              */
    const int64_t* off;
    const int* adj;         /* local column indices                         */
} BlockOp;

static void block_apply(const BlockOp* b, const double* x, double* y) {
    if (!b->two_stage) {
        for (int i = 0; i < b->k; i++) {
            double s = 0.0;
            for (int64_t p = b->off[i]; p < b->off[i + 1]; p++) s += x[b->adj[p]];
            y[i] = s;
        }
        return;
    }
    for (int i = 0; i < b->k; i++) y[i] = 0.0;
    for (int r = 0; r < b->nrows; r++) {
        double s = 0.0;
        for (int64_t p = b->off[r]; p < b->off[r + 1]; p++) s += x[b->adj[p]];
        for (int64_t p = b->off[r]; p < b->off[r + 1]; p++) y[b->adj[p]] += s;
    }
}

/* Perron vector of the k x k non-negative matrix M (row-major, destroyed):
 * M x = t x with t the eigenvalue of largest real part. */
static int dense_perron(double* M, int k, double* x) {
    if (k == 1) { x[0] = 1.0; return 1; }
    double* A = malloc((size_t)k * k * sizeof(double));
    double* wr = malloc((size_t)k * sizeof(double));
    double* wi = malloc((size_t)k * sizeof(double));
    double* VR = malloc((size_t)k * k * sizeof(double));
    int ok = 0;
    if (A && wr && wi && VR) {
        for (int i = 0; i < k; i++)                       /* to column-major */
            for (int j = 0; j < k; j++) A[(size_t)j * k + i] = M[(size_t)i * k + j];
        if (mat_lapack_dgeev(k, A, k, wr, wi, VR, k) == 0) {
            int best = -1;
            for (int j = 0; j < k; j++)
                if (fabs(wi[j]) <= 1e-12 * (fabs(wr[j]) + 1.0) && (best < 0 || wr[j] > wr[best]))
                    best = j;
            if (best >= 0) {
                double s = 0.0;
                for (int i = 0; i < k; i++) s += VR[(size_t)best * k + i];
                for (int i = 0; i < k; i++) {
                    double v = VR[(size_t)best * k + i] * (s < 0 ? -1.0 : 1.0);
                    x[i] = v > 0 ? v : 0.0;
                }
                ok = 1;
            }
        }
    }
    free(A); free(wr); free(wi); free(VR);
    return ok;
}

#define ARNOLDI_M 40
#define PERRON_TOL 1e-13

static int power_perron(const BlockOp* b, double* x) {
    int k = b->k;
    double* y = malloc((size_t)k * sizeof(double));
    if (!y) return 0;
    for (int i = 0; i < k; i++) x[i] = 1.0 / k;
    int ok = 0;
    for (int it = 0; it < 200000; it++) {
        block_apply(b, x, y);
        double t = 0.0, s = 0.0;
        for (int i = 0; i < k; i++) t += y[i];              /* sum(x) = 1 */
        double res = 0.0;
        for (int i = 0; i < k; i++) { double r = y[i] - t * x[i]; res += r * r; s += x[i] * x[i]; }
        if (t <= 0.0) break;
        if (sqrt(res) <= PERRON_TOL * t * sqrt(s)) { ok = 1; break; }
        for (int i = 0; i < k; i++) x[i] = (x[i] + y[i]) / (1.0 + t);   /* (I + M) shift */
    }
    free(y);
    return ok;
}

/* Perron vector of block b into x (non-negative, unnormalized). */
static int block_perron(const BlockOp* b, double* x) {
    int k = b->k;
    if (k == 1) { x[0] = 1.0; return 1; }
    if (k <= 64) {
        double* M = calloc((size_t)k * k, sizeof(double));
        double* e = calloc((size_t)k, sizeof(double));
        double* col = malloc((size_t)k * sizeof(double));
        int ok = M && e && col;
        for (int j = 0; j < k && ok; j++) {
            e[j] = 1.0;
            block_apply(b, e, col);
            e[j] = 0.0;
            for (int i = 0; i < k; i++) M[(size_t)i * k + j] = col[i];
        }
        if (ok) ok = dense_perron(M, k, x);
        free(M); free(e); free(col);
        if (ok) return 1;
        return power_perron(b, x);
    }
    int m = ARNOLDI_M < k ? ARNOLDI_M : k;
    double* V = malloc((size_t)(m + 1) * k * sizeof(double));
    double* H = malloc((size_t)(m + 1) * m * sizeof(double));   /* row-major (m+1) x m */
    double* Hc = malloc((size_t)m * m * sizeof(double));
    double* wr = malloc((size_t)m * sizeof(double));
    double* wi = malloc((size_t)m * sizeof(double));
    double* VR = malloc((size_t)m * m * sizeof(double));
    double* hv = malloc((size_t)(m + 1) * sizeof(double));
    if (!V || !H || !Hc || !wr || !wi || !VR || !hv) {
        free(V); free(H); free(Hc); free(wr); free(wi); free(VR); free(hv);
        return power_perron(b, x);
    }
    for (int i = 0; i < k; i++) x[i] = 1.0;
    int ok = 0, lapack_ok = 1;
    for (int restart = 0; restart < 2000 && !ok && lapack_ok; restart++) {
        double nx = 0.0;
        for (int i = 0; i < k; i++) nx += x[i] * x[i];
        nx = sqrt(nx);
        for (int i = 0; i < k; i++) V[i] = x[i] / nx;
        memset(H, 0, (size_t)(m + 1) * m * sizeof(double));
        int mm = m;
        for (int j = 0; j < m; j++) {
            double* w = V + (size_t)(j + 1) * k;
            block_apply(b, V + (size_t)j * k, w);
            double wn0 = 0.0;
            for (int i = 0; i < k; i++) wn0 += w[i] * w[i];
            for (int pass = 0; pass < 2; pass++) {         /* classical GS, twice */
#ifdef USE_LAPACK
                /* h = V_j w; w -= V_j^T h, as two BLAS matrix-vector products. */
                cblas_dgemv(CblasRowMajor, CblasNoTrans, j + 1, k, 1.0, V, k, w, 1, 0.0, hv, 1);
                cblas_dgemv(CblasRowMajor, CblasTrans, j + 1, k, -1.0, V, k, hv, 1, 1.0, w, 1);
                for (int q = 0; q <= j; q++) H[(size_t)q * m + j] += hv[q];
#else
                for (int q = 0; q <= j; q++) {
                    const double* vq = V + (size_t)q * k;
                    double h = 0.0;
                    for (int i = 0; i < k; i++) h += vq[i] * w[i];
                    H[(size_t)q * m + j] += h;
                    for (int i = 0; i < k; i++) w[i] -= h * vq[i];
                }
#endif
            }
            double wn = 0.0;
            for (int i = 0; i < k; i++) wn += w[i] * w[i];
            wn = sqrt(wn);
            H[(size_t)(j + 1) * m + j] = wn;
            if (wn <= 1e-14 * sqrt(wn0 > 0 ? wn0 : 1.0)) { mm = j + 1; break; }
            for (int i = 0; i < k; i++) w[i] /= wn;
        }
        for (int i = 0; i < mm; i++)
            for (int j = 0; j < mm; j++) Hc[(size_t)j * mm + i] = H[(size_t)i * m + j];
        if (mat_lapack_dgeev(mm, Hc, mm, wr, wi, VR, mm) != 0) { lapack_ok = 0; break; }
        int best = -1;
        for (int j = 0; j < mm; j++)
            if (fabs(wi[j]) <= 1e-10 * (fabs(wr[j]) + 1.0) && (best < 0 || wr[j] > wr[best]))
                best = j;
        if (best < 0) { lapack_ok = 0; break; }
        const double* y = VR + (size_t)best * mm;
        double s = 0.0;
        for (int i = 0; i < k; i++) {
            double v = 0.0;
            for (int q = 0; q < mm; q++) v += V[(size_t)q * k + i] * y[q];
            x[i] = v;
            s += v;
        }
        if (s < 0) for (int i = 0; i < k; i++) x[i] = -x[i];
        double theta = wr[best];
        double res = fabs(H[(size_t)mm * m + (mm - 1)] * y[mm - 1]);
        if (mm < m) res = 0.0;                          /* invariant subspace */
        if (theta > 0 && res <= PERRON_TOL * theta) ok = 1;
        if (theta <= 0) { lapack_ok = 0; break; }
    }
    free(V); free(H); free(Hc); free(wr); free(wi); free(VR); free(hv);
    if (ok) {
        for (int i = 0; i < k; i++) if (x[i] < 0) x[i] = 0.0;
        return 1;
    }
    return power_perron(b, x);
}

/* Combine per-block Perron vectors: block c (vertices with blk[v] == c) gets
 * its Perron vector scaled to total (size_c - 1) / sum (size - 1). rows_of(v)
 * yields, for one-stage blocks, the global row of v (columns are global
 * vertices, filtered to the block); for two-stage blocks the rows are the
 * global vertices u whose out-list lies in the block. Returns 1 and fills x
 * (length n), or 0 on failure. */
static int blocks_combine(int n, const int* blk, int nb, int two_stage,
                          const GmetCSR* rows, double* x) {
    size_t nn = (size_t)(n > 0 ? n : 1);
    int* size = calloc((size_t)(nb > 0 ? nb : 1), sizeof(int));
    int* start = calloc((size_t)nb + 1, sizeof(int));
    int* order = malloc(nn * sizeof(int));
    int* loc = malloc(nn * sizeof(int));
    if (!size || !start || !order || !loc) { free(size); free(start); free(order); free(loc); return 0; }
    for (int v = 0; v < n; v++) size[blk[v]]++;
    for (int c = 0; c < nb; c++) start[c + 1] = start[c] + size[c];
    int* fillp = malloc((size_t)(nb > 0 ? nb : 1) * sizeof(int));
    if (!fillp) { free(size); free(start); free(order); free(loc); return 0; }
    for (int c = 0; c < nb; c++) fillp[c] = start[c];
    for (int v = 0; v < n; v++) { loc[v] = fillp[blk[v]] - start[blk[v]]; order[fillp[blk[v]]++] = v; }
    free(fillp);
    /* Two-stage row owners: u belongs to the block of its first out-neighbour. */
    int* rstart = NULL; int* rorder = NULL;
    if (two_stage) {
        rstart = calloc((size_t)nb + 1, sizeof(int));
        rorder = malloc(nn * sizeof(int));
        int* rf = malloc((size_t)(nb > 0 ? nb : 1) * sizeof(int));
        if (!rstart || !rorder || !rf) { free(rstart); free(rorder); free(rf); free(size); free(start); free(order); free(loc); return 0; }
        for (int u = 0; u < n; u++)
            if (rows->off[u + 1] > rows->off[u]) rstart[blk[rows->adj[rows->off[u]]] + 1]++;
        for (int c = 0; c < nb; c++) rstart[c + 1] += rstart[c];
        for (int c = 0; c < nb; c++) rf[c] = rstart[c];
        for (int u = 0; u < n; u++)
            if (rows->off[u + 1] > rows->off[u]) rorder[rf[blk[rows->adj[rows->off[u]]]]++] = u;
        free(rf);
    }
    int64_t denom = 0;
    for (int c = 0; c < nb; c++) if (size[c] >= 2) denom += size[c] - 1;
    for (int v = 0; v < n; v++) x[v] = 0.0;
    int ok = 1;
    int64_t* boff = malloc((nn + 1) * sizeof(int64_t));
    int* badj = malloc((size_t)(rows->narcs + 1) * sizeof(int));
    double* bx = malloc(nn * sizeof(double));
    if (!boff || !badj || !bx) ok = 0;
    for (int c = 0; c < nb && ok; c++) {
        int k = size[c];
        if (k < 2) continue;
        BlockOp b;
        b.k = k; b.two_stage = two_stage; b.off = boff; b.adj = badj;
        int64_t q = 0;
        boff[0] = 0;
        if (!two_stage) {
            b.nrows = k;
            for (int i = 0; i < k; i++) {
                int v = order[start[c] + i];
                for (int64_t p = rows->off[v]; p < rows->off[v + 1]; p++) {
                    int u = rows->adj[p];
                    if (blk[u] == c) badj[q++] = loc[u];
                }
                boff[i + 1] = q;
            }
        } else {
            b.nrows = rstart[c + 1] - rstart[c];
            for (int r = 0; r < b.nrows; r++) {
                int u = rorder[rstart[c] + r];
                for (int64_t p = rows->off[u]; p < rows->off[u + 1]; p++) badj[q++] = loc[rows->adj[p]];
                boff[r + 1] = q;
            }
        }
        if (!block_perron(&b, bx)) { ok = 0; break; }
        double s = 0.0;
        for (int i = 0; i < k; i++) s += bx[i];
        double scale = s > 0 ? (double)(k - 1) / ((double)denom * s) : 0.0;
        for (int i = 0; i < k; i++) x[order[start[c] + i]] = bx[i] * scale;
    }
    free(boff); free(badj); free(bx);
    free(size); free(start); free(order); free(loc); free(rstart); free(rorder);
    return ok;
}

/* ---- EigenvectorCentrality ------------------------------------------------- */

Expr* builtin_eigenvector_centrality(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    int use_in = 1;
    if (argc == 2) {
        const Expr* m = res->data.function.args[1];
        if (m->type != EXPR_STRING) return NULL;
        if (strcmp(m->data.string, "In") == 0) use_in = 1;
        else if (strcmp(m->data.string, "Out") == 0) use_in = 0;
        else return NULL;
    }
    if (gmet_graph_kind(g) < 0) return NULL;
    Expr* cached = gmet_cache_get("EigenvectorCentrality", res);
    if (cached) return cached;
    GmetCSR out, pull;
    if (!gmet_csr_build(g, GMET_OUT, NULL, 0, &out)) return NULL;
    /* "In": x_v = sum over arcs u->v, so rows are the reversed arcs. */
    if (!gmet_csr_build(g, use_in ? GMET_IN : GMET_OUT, NULL, 0, &pull)) {
        gmet_csr_free(&out); return NULL;
    }
    int n = out.n;
    size_t nn = (size_t)(n > 0 ? n : 1);
    int* comp = malloc(nn * sizeof(int));
    double* x = malloc(nn * sizeof(double));
    Expr* result = NULL;
    if (comp && x) {
        int nc = gmet_scc(&out, comp);
        if (nc >= 0 && blocks_combine(n, comp, nc, 0, &pull, x))
            result = gmet_real_vector(x, n);
    }
    free(comp); free(x);
    gmet_csr_free(&out); gmet_csr_free(&pull);
    if (result) gmet_cache_put("EigenvectorCentrality", res, result);
    return result;
}

/* ---- KatzCentrality -------------------------------------------------------- */

Expr* builtin_katz_centrality(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 2 && argc != 3) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    int n = (int)g->data.function.args[0]->data.function.arg_count;
    size_t ne = g->data.function.args[1]->data.function.arg_count;
    double alpha;
    int azero;
    if (!read_real(res->data.function.args[1], &alpha, &azero)) return NULL;
    const Expr* bexpr = argc == 3 ? res->data.function.args[2] : NULL;
    int blist = bexpr && graph_is_list(bexpr);
    if (blist && (int)bexpr->data.function.arg_count != n) return NULL;
    size_t nn = (size_t)(n > 0 ? n : 1);
    double* b = malloc(nn * sizeof(double));
    if (!b) return NULL;
    for (int i = 0; i < n; i++) {
        const Expr* bi = !bexpr ? NULL : blist ? bexpr->data.function.args[i] : bexpr;
        double v = 1.0;
        if (bi && !read_real(bi, &v, NULL)) { free(b); return NULL; }
        b[i] = v;
    }
    if (ne == 0 || azero) {
        /* No propagation: Wolfram returns b itself, exact entries included. */
        free(b);
        Expr** it = malloc(nn * sizeof(Expr*));
        for (int i = 0; i < n; i++) {
            const Expr* bi = !bexpr ? NULL : blist ? bexpr->data.function.args[i] : bexpr;
            it[i] = bi ? expr_copy((Expr*)bi) : expr_new_integer(1);
        }
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), it, (size_t)n);
        free(it);
        return out;
    }
    Expr* cached = gmet_cache_get("KatzCentrality", res);
    if (cached) { free(b); return cached; }

    GmetCSR in;
    if (!gmet_csr_build(g, GMET_IN, NULL, 0, &in)) { free(b); return NULL; }
    double* x = malloc(nn * sizeof(double));
    double* xn = malloc(nn * sizeof(double));
    int converged = 0;
    if (x && xn) {
        memcpy(x, b, (size_t)n * sizeof(double));
        double prev = INFINITY;
        int grow = 0;
        for (int it = 0; it < 20000; it++) {
            double diff = 0.0, norm = 0.0;
            for (int v = 0; v < n; v++) {
                double s = 0.0;
                for (int64_t p = in.off[v]; p < in.off[v + 1]; p++) s += x[in.adj[p]];
                xn[v] = b[v] + alpha * s;
                diff += fabs(xn[v] - x[v]);
                norm += fabs(xn[v]);
            }
            double* t = x; x = xn; xn = t;
            if (!(norm < 1e300)) break;
            if (diff <= 1e-15 * norm) { converged = 1; break; }
            if (diff > prev) { if (++grow > 50) break; } else grow = 0;
            prev = diff;
        }
    }
    Expr* out = NULL;
    if (converged) out = gmet_real_vector(x, n);
    else if (x && xn && n <= 4000) {
        /* Dense solve of (I - alpha A^T) x = b, column-major. */
        double* M = calloc((size_t)n * n, sizeof(double));
        int* piv = malloc((size_t)n * sizeof(int));
        if (M && piv) {
            for (int v = 0; v < n; v++) {
                M[(size_t)v * n + v] += 1.0;
                for (int64_t p = in.off[v]; p < in.off[v + 1]; p++)
                    M[(size_t)in.adj[p] * n + v] -= alpha;     /* row v, column u */
            }
            memcpy(x, b, (size_t)n * sizeof(double));
            int info = mat_lapack_dgesv(n, 1, M, n, piv, x, n);
            int finite = info == 0;
            for (int v = 0; v < n && finite; v++) if (!isfinite(x[v])) finite = 0;
            if (finite) out = gmet_real_vector(x, n);
        }
        free(M); free(piv);
    }
    free(x); free(xn); free(b);
    gmet_csr_free(&in);
    if (out) gmet_cache_put("KatzCentrality", res, out);
    return out;
}

/* ---- HITSCentrality -------------------------------------------------------- */

static int uf_find(int* par, int x) {
    while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; }
    return x;
}

/* A^T A is block diagonal over the classes of "share an in-neighbour"; each
 * class is a block for blocks_combine (two-stage operator: the rows of B are
 * the out-lists of the class's in-neighbours). */
Expr* builtin_hits_centrality(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (gmet_graph_kind(g) < 0) return NULL;
    if (g->data.function.args[0]->data.function.arg_count == 0)       /* Wolfram: {} */
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    Expr* cached = gmet_cache_get("HITSCentrality", res);
    if (cached) return cached;
    GmetCSR out;
    if (!gmet_csr_build(g, GMET_OUT, NULL, 0, &out)) return NULL;
    int n = out.n;
    size_t nn = (size_t)(n > 0 ? n : 1);
    int* par = malloc(nn * sizeof(int));
    int* cls = calloc(nn, sizeof(int));
    double* h = malloc(nn * sizeof(double));
    double* a = malloc(nn * sizeof(double));
    Expr* result = NULL;
    if (par && cls && h && a) {
        for (int v = 0; v < n; v++) par[v] = v;
        for (int x = 0; x < n; x++)
            for (int64_t p = out.off[x] + 1; p < out.off[x + 1]; p++) {
                int r1 = uf_find(par, out.adj[out.off[x]]), r2 = uf_find(par, out.adj[p]);
                if (r1 != r2) par[r1] = r2;
            }
        int nc = 0;
        for (int v = 0; v < n; v++) if (uf_find(par, v) == v) cls[v] = nc++;
        for (int v = 0; v < n; v++) cls[v] = cls[uf_find(par, v)];
        if (blocks_combine(n, cls, nc, 1, &out, h)) {
            for (int u = 0; u < n; u++) {
                double s = 0.0;
                for (int64_t p = out.off[u]; p < out.off[u + 1]; p++) s += h[out.adj[p]];
                a[u] = s;
            }
            Expr* pair[2] = { gmet_real_vector(h, n), gmet_real_vector(a, n) };
            result = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
        }
    }
    free(par); free(cls); free(h); free(a);
    gmet_csr_free(&out);
    if (result) gmet_cache_put("HITSCentrality", res, result);
    return result;
}

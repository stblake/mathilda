/* gmet_distance.c - distance-based graph measures.
 *
 *   GraphDistanceMatrix[g]      all-pairs distances; [g, d] keeps those <= d
 *   GraphDistance[g, s]         distances from s to every vertex (the new
 *                               single-source form; [g, s, t] is delegated to
 *                               the existing builtin_graph_distance)
 *   VertexEccentricity[g, v]    largest distance from v to a vertex it reaches
 *   GraphDiameter / GraphRadius / GraphCenter / GraphPeriphery
 *   MeanGraphDistance[g]
 *   GraphDensity[g]
 *   KirchhoffMatrix[g]
 *
 * Semantics follow Mathematica 15 (checked with wolframscript on small graphs):
 *   - Distances follow edge direction; an undirected edge is usable both ways.
 *   - Unweighted graphs give exact Integers and Infinity; a graph carrying
 *     EdgeWeight gives machine Reals (Wolfram converts even integer weights),
 *     except that GraphDistance[g, s] reports the source itself as exact 0, as
 *     Wolfram does. A symbolic, complex or negative weight leaves the call
 *     unevaluated (Wolfram also refuses symbolic weights; negative weights,
 *     which Wolfram would route to Bellman-Ford, are not supported here).
 *   - Unweighted: VertexEccentricity measures over the vertices v can reach,
 *     so it is finite on a disconnected graph, and GraphDiameter/GraphRadius/
 *     MeanGraphDistance are Infinity, GraphCenter/GraphPeriphery {}, unless
 *     the graph is strongly connected (connected, when undirected).
 *   - Weighted (Wolfram's other rule, verified): an unreachable vertex makes
 *     VertexEccentricity Infinity, and radius/center/periphery are taken over
 *     those eccentricities -- e.g. a weighted digraph that is not strongly
 *     connected can have a finite radius, and a disconnected weighted graph
 *     has every vertex in both its center and its periphery. Diameter and
 *     MeanGraphDistance are Infinity whenever some pair is unreachable.
 *   - MeanGraphDistance averages over ordered pairs of distinct vertices and is
 *     exact (Integer/Rational) when unweighted. Undefined for one vertex.
 *   - GraphDensity = (directed edges + 2 undirected edges) / (n (n - 1)).
 *   - KirchhoffMatrix = D - A with D the number of edges incident to each
 *     vertex and A the (directed) adjacency matrix; weights are ignored.
 *     Wolfram returns a SparseArray; Mathilda has none, so this is the dense
 *     (packed Integer) matrix -- Normal of Wolfram's answer.
 *
 * Algorithms: bit-parallel multi-source BFS (gmet_core.c) for unweighted
 * all-pairs work -- 256 sources per adjacency sweep, batches spread over the
 * thread team -- and binary-heap Dijkstra per source for weighted graphs.
 * Diameter/radius/center/periphery/mean distance reduce from the cached
 * per-source summary, after an O(n + m) strong-connectivity test that settles
 * the disconnected case without any all-pairs work.
 */

#include "graph_metrics.h"
#include "graph.h"
#include "expr.h"
#include "eval.h"
#include "pack.h"
#include "arithmetic.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---- helpers --------------------------------------------------------------- */

static Expr* list_of(Expr** items, size_t n) {
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, n);
    return out;
}

/* A distance cutoff argument: a non-negative real number, or Infinity. Returns
 * 1 and sets *d on success. */
static int read_cutoff(const Expr* e, double* d) {
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity) {
        *d = INFINITY; return 1;
    }
    double x = graph_weight_to_double(e);
    if (!(x >= 0.0)) return 0;
    *d = x;
    return 1;
}

/* True iff the OUT arcs of g make it strongly connected (n >= 1). */
static int strongly_connected(const GmetCSR* c) {
    if (c->n <= 1) return 1;
    int* comp = malloc((size_t)c->n * sizeof(int));
    if (!comp) return -1;
    int k = gmet_scc(c, comp);
    free(comp);
    return k < 0 ? -1 : k == 1;
}

/* ---- GraphDistanceMatrix ------------------------------------------------- */

typedef struct {
    const GmetCSR* g;
    int64_t* D;                /* unweighted: n*n, -1 = unreached            */
    double* Dw;                /* weighted: n*n                              */
    int maxlevel;
    GmetMSBFS* bfs[32];
    GmetHeap* heap[32];
} GdmCtx;

typedef struct { int64_t* row0; int n; int base; } GdmBatch;

static void gdm_visit(void* ctx, int level, int v, const uint64_t* lanes) {
    GdmBatch* b = (GdmBatch*)ctx;
    int64_t* row = b->row0 + (size_t)v * (size_t)b->n + b->base;
    for (int k = 0; k < GMET_MSBFS_WORDS; k++) {
        uint64_t x = lanes[k];
        while (x) {
            int i = (k << 6) + GMET_CTZ64(x);
            x &= x - 1;
            row[i] = level;
        }
    }
}

static void gdm_bfs_work(void* vctx, int tid, int64_t lo, int64_t hi) {
    GdmCtx* c = (GdmCtx*)vctx;
    int n = c->g->n;
    int src[GMET_MSBFS_LANES];
    for (int64_t blk = lo; blk < hi; blk++) {
        int base = (int)(blk * GMET_MSBFS_LANES);
        int ns = n - base < GMET_MSBFS_LANES ? n - base : GMET_MSBFS_LANES;
        for (int i = 0; i < ns; i++) src[i] = base + i;
        GdmBatch b = { c->D, n, base };
        gmet_msbfs_run(c->bfs[tid], c->g, src, ns, c->maxlevel, gdm_visit, &b);
    }
}

static void gdm_dijkstra_work(void* vctx, int tid, int64_t lo, int64_t hi) {
    GdmCtx* c = (GdmCtx*)vctx;
    size_t n = (size_t)c->g->n;
    for (int64_t s = lo; s < hi; s++)
        gmet_dijkstra(c->g, (int)s, c->Dw + (size_t)s * n, NULL, NULL, c->heap[tid]);
}

/* Unpacked matrix with Infinity entries (and exact/real entries otherwise). */
static Expr* matrix_from_int(const int64_t* D, int n) {
    Expr** rows = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    Expr* inf = expr_new_symbol(SYM_Infinity);
    for (int i = 0; i < n; i++) {
        Expr** row = malloc((size_t)n * sizeof(Expr*));
        for (int j = 0; j < n; j++) {
            int64_t d = D[(size_t)i * n + j];
            row[j] = d < 0 ? expr_copy(inf) : expr_new_integer(d);
        }
        Expr* r = pack_offer(list_of(row, (size_t)n));
        free(row);
        rows[i] = r;
    }
    expr_free(inf);
    Expr* out = list_of(rows, (size_t)n);
    free(rows);
    return out;
}

static Expr* matrix_from_real(const double* D, int n) {
    Expr** rows = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    Expr* inf = expr_new_symbol(SYM_Infinity);
    for (int i = 0; i < n; i++) {
        Expr** row = malloc((size_t)n * sizeof(Expr*));
        for (int j = 0; j < n; j++) {
            double d = D[(size_t)i * n + j];
            row[j] = isinf(d) ? expr_copy(inf) : expr_new_real(d);
        }
        Expr* r = pack_offer(list_of(row, (size_t)n));
        free(row);
        rows[i] = r;
    }
    expr_free(inf);
    Expr* out = list_of(rows, (size_t)n);
    free(rows);
    return out;
}

Expr* builtin_graph_distance_matrix(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    double cutoff = INFINITY;
    if (argc == 2 && !read_cutoff(res->data.function.args[1], &cutoff)) return NULL;
    int n = (int)g->data.function.args[0]->data.function.arg_count;
    if (n == 0) return NULL;                       /* Wolfram: unevaluated */

    Expr* cached = gmet_cache_get("GraphDistanceMatrix", res);
    if (cached) return cached;

    double* ew = NULL;
    int wk = gmet_edge_weights(g, &ew);
    if (wk < 0) return NULL;
    size_t nn = (size_t)n * (size_t)n;
    int64_t dims[2] = { n, n };
    GdmCtx c;
    memset(&c, 0, sizeof(c));
    Expr* out = NULL;

    if (wk == 0) {
        /* BFS over the REVERSED arcs from target t reaches w at level L exactly
         * when D[w][t] = L, so each batch writes contiguous runs of row w. */
        GmetCSR csr;
        if (!gmet_csr_build(g, GMET_IN, NULL, 0, &csr)) return NULL;
        void* buf = NULL;
        Expr* nd = ndbuild_open(2, dims, NDT_INT64, &buf);
        int64_t* D = nd ? (int64_t*)buf : malloc(nn * sizeof(int64_t));
        if (!D) { gmet_csr_free(&csr); return NULL; }
        for (size_t i = 0; i < nn; i++) D[i] = -1;
        for (int i = 0; i < n; i++) D[(size_t)i * n + i] = 0;
        c.g = &csr; c.D = D;
        c.maxlevel = isinf(cutoff) ? -1 : (int)floor(cutoff);
        if (c.maxlevel >= 0 || isinf(cutoff)) {
            int64_t nblk = (n + GMET_MSBFS_LANES - 1) / GMET_MSBFS_LANES;
            int nt = gmet_thread_count((double)nblk * (double)(csr.narcs + n) * 8.0);
            if (nt > nblk) nt = (int)nblk;
            nt = gmet_cap_threads(nt, (double)n * (3 * 8 * GMET_MSBFS_WORDS + 14));
            int ok = 1;
            for (int t = 0; t < nt; t++) if (!(c.bfs[t] = gmet_msbfs_new(n))) ok = 0;
            if (ok) gmet_parallel_for(nblk, 1, nt, gdm_bfs_work, &c);
            for (int t = 0; t < nt; t++) gmet_msbfs_free(c.bfs[t]);
            if (!ok) { gmet_csr_free(&csr); if (nd) expr_free(nd); else free(D); return NULL; }
        }
        gmet_csr_free(&csr);
        int complete = 1;
        for (size_t i = 0; i < nn && complete; i++) if (D[i] < 0) complete = 0;
        if (complete && nd) out = nd;
        else {
            out = matrix_from_int(D, n);
            if (nd) expr_free(nd); else free(D);
        }
    } else {
        GmetCSR csr;
        int okb = gmet_csr_build(g, GMET_OUT, ew, 0, &csr);
        free(ew);
        if (!okb) return NULL;
        void* buf = NULL;
        Expr* nd = ndbuild_open(2, dims, NDT_FLOAT64, &buf);
        double* Dw = nd ? (double*)buf : malloc(nn * sizeof(double));
        if (!Dw) { gmet_csr_free(&csr); return NULL; }
        c.g = &csr; c.Dw = Dw;
        int nt = gmet_cap_threads(gmet_thread_count((double)n * (double)(csr.narcs * 4 + n)),
                                  (double)(csr.narcs + n) * 12);
        int ok = 1;
        for (int t = 0; t < nt; t++) if (!(c.heap[t] = gmet_heap_new(csr.narcs + n + 1))) ok = 0;
        if (ok) gmet_parallel_for(n, 8, nt, gdm_dijkstra_work, &c);
        for (int t = 0; t < nt; t++) gmet_heap_free(c.heap[t]);
        gmet_csr_free(&csr);
        if (!ok) { if (nd) expr_free(nd); else free(Dw); return NULL; }
        int complete = 1;
        for (size_t i = 0; i < nn; i++) {
            if (Dw[i] > cutoff) Dw[i] = INFINITY;
            if (isinf(Dw[i])) complete = 0;
        }
        if (complete && nd) out = nd;
        else {
            out = matrix_from_real(Dw, n);
            if (nd) expr_free(nd); else free(Dw);
        }
    }
    gmet_cache_put("GraphDistanceMatrix", res, out);
    return out;
}

/* ---- GraphDistance[g, s] -------------------------------------------------- */

/* Single-source distances from vertex index s. Unweighted: BFS, int64 with -1
 * for unreached; weighted: Dijkstra into dw. */
static int single_source(const Expr* g, int s, int wk, const double* ew,
                         int64_t* di, double* dw) {
    GmetCSR csr;
    if (!gmet_csr_build(g, GMET_OUT, wk ? ew : NULL, 0, &csr)) return 0;
    int n = csr.n;
    if (!wk) {
        int* q = malloc((size_t)n * sizeof(int));
        if (!q) { gmet_csr_free(&csr); return 0; }
        for (int i = 0; i < n; i++) di[i] = -1;
        int h = 0, t = 0;
        di[s] = 0; q[t++] = s;
        while (h < t) {
            int u = q[h++];
            for (int64_t p = csr.off[u]; p < csr.off[u + 1]; p++) {
                int v = csr.adj[p];
                if (di[v] < 0) { di[v] = di[u] + 1; q[t++] = v; }
            }
        }
        free(q);
    } else {
        GmetHeap* hp = gmet_heap_new(csr.narcs + n + 1);
        if (!hp) { gmet_csr_free(&csr); return 0; }
        gmet_dijkstra(&csr, s, dw, NULL, NULL, hp);
        gmet_heap_free(hp);
    }
    gmet_csr_free(&csr);
    return 1;
}

Expr* builtin_gmet_graph_distance(Expr* res) {
    if (res->data.function.arg_count != 2) return builtin_graph_distance(res);
    const Expr* g = res->data.function.args[0];
    int s = graph_vertex_position(g, res->data.function.args[1]);
    if (s < 0) return NULL;
    int n = (int)g->data.function.args[0]->data.function.arg_count;
    double* ew = NULL;
    int wk = gmet_edge_weights(g, &ew);
    if (wk < 0) return NULL;
    int64_t* di = malloc((size_t)n * sizeof(int64_t));
    double* dw = malloc((size_t)n * sizeof(double));
    if (!di || !dw || !single_source(g, s, wk, ew, di, dw)) {
        free(di); free(dw); free(ew); return NULL;
    }
    free(ew);
    Expr* out;
    int complete = 1;
    if (!wk) {
        for (int i = 0; i < n; i++) if (di[i] < 0) complete = 0;
        if (complete) out = gmet_int_vector(di, n);
        else {
            Expr** it = malloc((size_t)n * sizeof(Expr*));
            for (int i = 0; i < n; i++)
                it[i] = di[i] < 0 ? expr_new_symbol(SYM_Infinity) : expr_new_integer(di[i]);
            out = list_of(it, (size_t)n);
            free(it);
        }
    } else {
        /* Wolfram: machine reals, but the source itself is an exact 0. */
        Expr** it = malloc((size_t)n * sizeof(Expr*));
        for (int i = 0; i < n; i++)
            it[i] = i == s ? expr_new_integer(0) : gmet_distance_value(dw[i], 1);
        out = list_of(it, (size_t)n);
        free(it);
    }
    free(di); free(dw);
    return out;
}

/* ---- VertexEccentricity ---------------------------------------------------- */

Expr* builtin_vertex_eccentricity(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    int s = graph_vertex_position(g, res->data.function.args[1]);
    if (s < 0) return NULL;
    int n = (int)g->data.function.args[0]->data.function.arg_count;
    double* ew = NULL;
    int wk = gmet_edge_weights(g, &ew);
    if (wk < 0) return NULL;
    int64_t* di = malloc((size_t)n * sizeof(int64_t));
    double* dw = malloc((size_t)n * sizeof(double));
    if (!di || !dw || !single_source(g, s, wk, ew, di, dw)) {
        free(di); free(dw); free(ew); return NULL;
    }
    free(ew);
    Expr* out;
    if (!wk) {
        int64_t e = 0;
        for (int i = 0; i < n; i++) if (di[i] > e) e = di[i];
        out = expr_new_integer(e);
    } else {
        /* Weighted: Wolfram measures over ALL vertices, so an unreachable one
         * makes the eccentricity Infinity (unweighted: reachable ones only). */
        double e = 0.0;
        for (int i = 0; i < n; i++) if (dw[i] > e) e = dw[i];
        out = gmet_distance_value(e, 1);
    }
    free(di); free(dw);
    return out;
}

/* ---- Diameter / radius / center / periphery / mean distance --------------- */

enum { EX_DIAMETER, EX_RADIUS, EX_CENTER, EX_PERIPHERY, EX_MEAN };

/* Relative tolerance for "equal eccentricity" on weighted graphs, so that two
 * sums of the same weights taken in different orders compare equal. */
#define ECC_TIE_TOL 1e-12

static const char* const EX_NAMES[] = {
    "GraphDiameter", "GraphRadius", "GraphCenter", "GraphPeriphery", "MeanGraphDistance"
};

static Expr* extremal_compute(Expr* res, int what);

static Expr* extremal(Expr* res, int what) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* cached = gmet_cache_get(EX_NAMES[what], res);
    if (cached) return cached;
    Expr* out = extremal_compute(res, what);
    if (out) gmet_cache_put(EX_NAMES[what], res, out);
    return out;
}

static Expr* extremal_compute(Expr* res, int what) {
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    int n = (int)g->data.function.args[0]->data.function.arg_count;
    double* ew = NULL;
    int wk = gmet_edge_weights(g, &ew);
    free(ew);
    if (wk < 0) return NULL;

    if (n == 0) {
        if (what == EX_DIAMETER || what == EX_RADIUS) return expr_new_integer(0);
        if (what == EX_MEAN) return NULL;
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }
    if (what == EX_MEAN && n == 1) return NULL;   /* Wolfram: unevaluated */

    GmetCSR csr;
    if (!gmet_csr_build(g, GMET_OUT, NULL, 0, &csr)) return NULL;
    int sc = strongly_connected(&csr);
    gmet_csr_free(&csr);
    if (sc < 0) return NULL;
    if (!sc && what == EX_MEAN) return expr_new_symbol(SYM_Infinity);
    if (!sc && !wk) {
        if (what == EX_CENTER || what == EX_PERIPHERY)
            return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        return expr_new_symbol(SYM_Infinity);
    }

    int status;
    const GmetDistSummary* S = gmet_dist_summary(g, &status);
    if (!S) return NULL;
    /* Weighted graphs use Wolfram's weighted eccentricity, Infinity for a
     * vertex that does not reach every other one (see VertexEccentricity). */
    double* ecc = malloc((size_t)n * sizeof(double));
    if (!ecc) return NULL;
    double lo = INFINITY, hi = 0.0, tot = 0.0;
    for (int i = 0; i < n; i++) {
        ecc[i] = (S->weighted && S->reach[i] < n - 1) ? INFINITY : S->ecc[i];
        if (ecc[i] < lo) lo = ecc[i];
        if (ecc[i] > hi) hi = ecc[i];
        tot += S->sum[i];
    }
    Expr* out = NULL;
    switch (what) {
    case EX_DIAMETER: out = gmet_distance_value(hi, S->weighted); break;
    case EX_RADIUS:   out = gmet_distance_value(lo, S->weighted); break;
    case EX_MEAN:
        if (S->weighted) out = expr_new_real(tot / ((double)n * (double)(n - 1)));
        else out = make_rational((int64_t)tot, (int64_t)n * (int64_t)(n - 1));
        break;
    default: {
        double target = what == EX_CENTER ? lo : hi;
        double tol = (S->weighted && !isinf(target)) ? ECC_TIE_TOL * (fabs(target) + 1.0) : 0.0;
        unsigned char* flag = calloc((size_t)n, 1);
        if (flag) {
            for (int i = 0; i < n; i++)
                flag[i] = isinf(target) ? isinf(ecc[i]) : fabs(ecc[i] - target) <= tol;
            out = gmet_vertex_subset(g, flag);
            free(flag);
        }
    }
    }
    free(ecc);
    return out;
}

Expr* builtin_graph_diameter(Expr* res)      { return extremal(res, EX_DIAMETER); }
Expr* builtin_graph_radius(Expr* res)        { return extremal(res, EX_RADIUS); }
Expr* builtin_graph_center(Expr* res)        { return extremal(res, EX_CENTER); }
Expr* builtin_graph_periphery(Expr* res)     { return extremal(res, EX_PERIPHERY); }
Expr* builtin_mean_graph_distance(Expr* res) { return extremal(res, EX_MEAN); }

/* ---- GraphDensity --------------------------------------------------------- */

Expr* builtin_graph_density(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    long nd = graph_directed_edge_count(g);
    if (nd < 0) return NULL;
    int64_t n = (int64_t)g->data.function.args[0]->data.function.arg_count;
    int64_t ne = (int64_t)g->data.function.args[1]->data.function.arg_count;
    if (n < 2) return NULL;                         /* Wolfram: unevaluated */
    return make_rational((int64_t)nd + 2 * (ne - (int64_t)nd), n * (n - 1));
}

/* ---- KirchhoffMatrix ------------------------------------------------------ */

Expr* builtin_kirchhoff_matrix(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    const int *eu, *ev;
    const unsigned char* edir;
    if (!graph_edge_indices(g, &eu, &ev, &edir)) return NULL;
    int n = (int)g->data.function.args[0]->data.function.arg_count;
    size_t ne = g->data.function.args[1]->data.function.arg_count;
    if (n == 0) return NULL;
    size_t nn = (size_t)n * (size_t)n;
    int64_t dims[2] = { n, n };
    void* buf = NULL;
    Expr* nd = ndbuild_open(2, dims, NDT_INT64, &buf);
    int64_t* K = nd ? (int64_t*)buf : calloc(nn, sizeof(int64_t));
    if (!K) return NULL;
    if (nd) memset(K, 0, nn * sizeof(int64_t));
    for (size_t k = 0; k < ne; k++) {
        size_t a = (size_t)eu[k], b = (size_t)ev[k];
        K[a * n + a]++; K[b * n + b]++;
        K[a * n + b] = -1;
        if (!edir[k]) K[b * n + a] = -1;
    }
    if (nd) return nd;
    Expr** rows = malloc((size_t)n * sizeof(Expr*));
    for (int i = 0; i < n; i++) {
        Expr** row = malloc((size_t)n * sizeof(Expr*));
        for (int j = 0; j < n; j++) row[j] = expr_new_integer(K[(size_t)i * n + j]);
        rows[i] = list_of(row, (size_t)n);
        free(row);
    }
    free(K);
    Expr* out = list_of(rows, (size_t)n);
    free(rows);
    return out;
}

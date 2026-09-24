/* gmet_centrality.c - path-based centralities.
 *
 *   DegreeCentrality[g] / [g, "In"] / [g, "Out"]
 *   ClosenessCentrality[g]
 *   EccentricityCentrality[g]
 *   BetweennessCentrality[g]
 *   EdgeBetweennessCentrality[g]
 *
 * Semantics follow Mathematica 15 (verified with wolframscript):
 *   - DegreeCentrality gives exact Integers. On an undirected graph every mode
 *     is the vertex degree. Otherwise an undirected edge counts as one arc in
 *     each direction, "In"/"Out" count incoming/outgoing arcs and the default
 *     is their sum -- so in a MIXED graph an undirected edge adds 2 to the
 *     total of each endpoint, as Wolfram does. Weights are ignored.
 *   - ClosenessCentrality[g][[v]] = r / (sum of distances from v to the r
 *     vertices it reaches), 0 when it reaches none. Weighted (EdgeWeight used
 *     as lengths). Machine reals.
 *   - EccentricityCentrality[g][[v]] = 1 / VertexEccentricity[g, v], 0 when the
 *     eccentricity is 0. Weighted. Machine reals.
 *   - BetweennessCentrality is the unnormalized Freeman betweenness; weights
 *     are IGNORED (Wolfram gives the same answer with or without EdgeWeight).
 *     On an undirected graph each unordered pair is counted once; on a
 *     directed graph each ordered pair. Wolfram's value on MIXED graphs does
 *     not match any standard definition (e.g. {1<->2, 2->3} gives vertex 2 a
 *     betweenness of 0), so a mixed graph is left unevaluated rather than
 *     guessed at.
 *   - EdgeBetweennessCentrality sums over ORDERED pairs even on undirected
 *     graphs (Wolfram: 2 for the only edge of a K2), in EdgeList order, and
 *     DOES use edge weights as lengths (ties within a relative 1e-12 share
 *     paths). Mixed graphs follow the arc semantics (undirected = both ways).
 *
 * Algorithms: closeness/eccentricity reduce from the cached multi-source-BFS /
 * Dijkstra distance summary (gmet_core.c). Betweenness is Brandes' algorithm
 * (O(nm) unweighted, O(nm + n^2 log n) weighted), sources spread over the
 * thread team with per-thread accumulators merged at the end, successor-
 * scanning accumulation (no predecessor lists).
 */

#include "graph_metrics.h"
#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---- DegreeCentrality ------------------------------------------------------ */

Expr* builtin_degree_centrality(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    int mode = 0;                                   /* 0 all, 1 in, 2 out */
    if (argc == 2) {
        const Expr* m = res->data.function.args[1];
        if (m->type != EXPR_STRING) return NULL;
        if (strcmp(m->data.string, "In") == 0) mode = 1;
        else if (strcmp(m->data.string, "Out") == 0) mode = 2;
        else return NULL;
    }
    int kind = gmet_graph_kind(g);
    if (kind < 0) return NULL;
    const int *eu, *ev;
    const unsigned char* edir;
    graph_edge_indices(g, &eu, &ev, &edir);
    int n = (int)g->data.function.args[0]->data.function.arg_count;
    size_t ne = g->data.function.args[1]->data.function.arg_count;
    int64_t* in = calloc((size_t)(n > 0 ? n : 1), sizeof(int64_t));
    int64_t* out = calloc((size_t)(n > 0 ? n : 1), sizeof(int64_t));
    if (!in || !out) { free(in); free(out); return NULL; }
    for (size_t k = 0; k < ne; k++) {
        out[eu[k]]++; in[ev[k]]++;
        if (!edir[k]) { out[ev[k]]++; in[eu[k]]++; }
    }
    int64_t* d = mode == 1 ? in : out;
    if (mode == 0) {
        if (kind == GMET_KIND_UNDIRECTED) d = out;   /* out == in == degree */
        else for (int i = 0; i < n; i++) out[i] += in[i];
    }
    Expr* r = gmet_int_vector(d, n);
    free(in); free(out);
    return r;
}

/* ---- Closeness & eccentricity centrality ----------------------------------- */

static Expr* summary_centrality(Expr* res, int ecc) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    int status;
    const GmetDistSummary* S = gmet_dist_summary(g, &status);
    if (!S) return NULL;
    int n = S->n;
    double* v = malloc((size_t)(n > 0 ? n : 1) * sizeof(double));
    if (!v) return NULL;
    for (int i = 0; i < n; i++) {
        if (ecc) v[i] = S->ecc[i] > 0 ? 1.0 / S->ecc[i] : 0.0;
        else v[i] = (S->reach[i] > 0 && S->sum[i] > 0) ? (double)S->reach[i] / S->sum[i] : 0.0;
    }
    Expr* out = gmet_real_vector(v, n);
    free(v);
    return out;
}

Expr* builtin_closeness_centrality(Expr* res)    { return summary_centrality(res, 0); }
Expr* builtin_eccentricity_centrality(Expr* res) { return summary_centrality(res, 1); }

/* ---- Brandes -------------------------------------------------------------- */

typedef struct {
    const GmetCSR* g;
    int edges;            /* accumulate per edge (else per vertex)            */
    int nacc;             /* accumulator length                               */
    double* acc[32];      /* per-thread accumulators                          */
    /* per-thread scratch */
    int* dist[32];        /* BFS                                              */
    double* dw[32];       /* Dijkstra                                         */
    double* sigma[32];
    double* delta[32];
    int* order[32];
    int* pos[32];
    GmetHeap* heap[32];
} BrandesCtx;

static void brandes_bfs_work(void* vctx, int tid, int64_t lo, int64_t hi) {
    BrandesCtx* c = (BrandesCtx*)vctx;
    const GmetCSR* g = c->g;
    int n = g->n;
    int* dist = c->dist[tid];
    double* sigma = c->sigma[tid];
    double* delta = c->delta[tid];
    int* S = c->order[tid];
    double* acc = c->acc[tid];
    const int64_t* off = g->off;
    const int* adj = g->adj;
    for (int i = 0; i < n; i++) dist[i] = -1;
    for (int64_t s = lo; s < hi; s++) {
        int ns = 0, head = 0;
        dist[s] = 0; sigma[s] = 1.0; S[ns++] = (int)s;
        while (head < ns) {
            int v = S[head++];
            int dv1 = dist[v] + 1;
            double sv = sigma[v];
            for (int64_t p = off[v]; p < off[v + 1]; p++) {
                int w = adj[p];
                if (dist[w] < 0) { dist[w] = dv1; sigma[w] = sv; S[ns++] = w; }
                else if (dist[w] == dv1) sigma[w] += sv;
            }
        }
        for (int k = ns - 1; k >= 0; k--) {
            int v = S[k];
            int dv1 = dist[v] + 1;
            double coeff = sigma[v];
            double dsum = 0.0;
            for (int64_t p = off[v]; p < off[v + 1]; p++) {
                int w = adj[p];
                if (dist[w] == dv1) {
                    double cr = coeff / sigma[w] * (1.0 + delta[w]);
                    dsum += cr;
                    if (c->edges) acc[g->eid[p]] += cr;
                }
            }
            delta[v] = dsum;
            if (!c->edges && v != s) acc[v] += dsum;
        }
        for (int k = 0; k < ns; k++) { dist[S[k]] = -1; delta[S[k]] = 0.0; }
    }
}

#define BRANDES_TIE 1e-12

static void brandes_dijkstra_work(void* vctx, int tid, int64_t lo, int64_t hi) {
    BrandesCtx* c = (BrandesCtx*)vctx;
    const GmetCSR* g = c->g;
    int n = g->n;
    double* dist = c->dw[tid];
    double* sigma = c->sigma[tid];
    double* delta = c->delta[tid];
    int* S = c->order[tid];
    int* pos = c->pos[tid];
    double* acc = c->acc[tid];
    GmetHeap* h = c->heap[tid];
    for (int64_t s = lo; s < hi; s++) {
        /* Dijkstra with path counting: sigma[v] accumulates from every vertex
         * settled before v that lies on a shortest path to it. */
        int ns;
        for (int i = 0; i < n; i++) { pos[i] = -1; sigma[i] = 0.0; delta[i] = 0.0; }
        gmet_dijkstra(g, (int)s, dist, S, &ns, h);
        for (int k = 0; k < ns; k++) pos[S[k]] = k;
        sigma[s] = 1.0;
        for (int k = 0; k < ns; k++) {
            int v = S[k];
            double dv = dist[v];
            for (int64_t p = g->off[v]; p < g->off[v + 1]; p++) {
                int w = g->adj[p];
                if (pos[w] <= k) continue;               /* settled no later than v */
                double nd = dv + g->w[p];
                if (fabs(nd - dist[w]) <= BRANDES_TIE * (dist[w] > 1.0 ? dist[w] : 1.0))
                    sigma[w] += sigma[v];
            }
        }
        for (int k = ns - 1; k >= 0; k--) {
            int v = S[k];
            double dv = dist[v];
            double dsum = 0.0;
            for (int64_t p = g->off[v]; p < g->off[v + 1]; p++) {
                int w = g->adj[p];
                if (pos[w] <= k) continue;
                double nd = dv + g->w[p];
                if (fabs(nd - dist[w]) <= BRANDES_TIE * (dist[w] > 1.0 ? dist[w] : 1.0)) {
                    double cr = sigma[v] / sigma[w] * (1.0 + delta[w]);
                    dsum += cr;
                    if (c->edges) acc[g->eid[p]] += cr;
                }
            }
            delta[v] = dsum;
            if (!c->edges && v != s) acc[v] += dsum;
        }
    }
}

/* Runs Brandes over every source; returns a malloc'd accumulator of length
 * c->nacc (vertices or edges), or NULL on allocation failure. */
static double* brandes_run(const GmetCSR* g, int edges, int nacc) {
    int n = g->n;
    BrandesCtx c;
    memset(&c, 0, sizeof(c));
    c.g = g; c.edges = edges; c.nacc = nacc;
    double work = (double)n * (double)(g->narcs * (g->w ? 6 : 2) + n);
    int nt = gmet_cap_threads(gmet_thread_count(work),
                              8.0 * nacc + 32.0 * n + (g->w ? 12.0 * (double)(g->narcs + n) : 0.0));
    size_t nn = (size_t)(n > 0 ? n : 1), na = (size_t)(nacc > 0 ? nacc : 1);
    int ok = 1;
    for (int t = 0; t < nt; t++) {
        c.acc[t] = calloc(na, sizeof(double));
        c.sigma[t] = malloc(nn * sizeof(double));
        c.delta[t] = calloc(nn, sizeof(double));
        c.order[t] = malloc(nn * sizeof(int));
        if (g->w) {
            c.dw[t] = malloc(nn * sizeof(double));
            c.pos[t] = malloc(nn * sizeof(int));
            c.heap[t] = gmet_heap_new(g->narcs + n + 1);
            if (!c.dw[t] || !c.pos[t] || !c.heap[t]) ok = 0;
        } else {
            c.dist[t] = malloc(nn * sizeof(int));
            if (!c.dist[t]) ok = 0;
        }
        if (!c.acc[t] || !c.sigma[t] || !c.delta[t] || !c.order[t]) ok = 0;
    }
    if (ok) {
        int64_t chunk = n / (nt * 16) + 1;
        if (chunk > 64) chunk = 64;
        gmet_parallel_for(n, chunk, nt, g->w ? brandes_dijkstra_work : brandes_bfs_work, &c);
        for (int t = 1; t < nt; t++)
            for (int i = 0; i < nacc; i++) c.acc[0][i] += c.acc[t][i];
    }
    double* out = ok ? c.acc[0] : NULL;
    for (int t = 0; t < nt; t++) {
        if (t > 0 || !ok) free(c.acc[t]);
        free(c.sigma[t]); free(c.delta[t]); free(c.order[t]);
        free(c.dist[t]); free(c.dw[t]); free(c.pos[t]); gmet_heap_free(c.heap[t]);
    }
    return out;
}

Expr* builtin_betweenness_centrality(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    int kind = gmet_graph_kind(g);
    if (kind < 0 || kind == GMET_KIND_MIXED) return NULL;
    Expr* cached = gmet_cache_get("BetweennessCentrality", res);
    if (cached) return cached;
    GmetCSR csr;
    if (!gmet_csr_build(g, GMET_OUT, NULL, 0, &csr)) return NULL;
    int n = csr.n;
    double* cb = brandes_run(&csr, 0, n);
    gmet_csr_free(&csr);
    if (!cb) return NULL;
    if (kind == GMET_KIND_UNDIRECTED) for (int i = 0; i < n; i++) cb[i] *= 0.5;
    Expr* out = gmet_real_vector(cb, n);
    free(cb);
    gmet_cache_put("BetweennessCentrality", res, out);
    return out;
}

Expr* builtin_edge_betweenness_centrality(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    Expr* cached = gmet_cache_get("EdgeBetweennessCentrality", res);
    if (cached) return cached;
    double* ew = NULL;
    int wk = gmet_edge_weights(g, &ew);
    if (wk < 0) return NULL;
    GmetCSR csr;
    int okb = gmet_csr_build(g, GMET_OUT, ew, 1, &csr);
    free(ew);
    if (!okb) return NULL;
    int ne = (int)g->data.function.args[1]->data.function.arg_count;
    double* eb = brandes_run(&csr, 1, ne);
    gmet_csr_free(&csr);
    if (!eb) return NULL;
    Expr* out = gmet_real_vector(eb, ne);
    free(eb);
    gmet_cache_put("EdgeBetweennessCentrality", res, out);
    return out;
}

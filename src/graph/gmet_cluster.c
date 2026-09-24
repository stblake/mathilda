/* gmet_cluster.c - triangle counting and clustering coefficients.
 *
 *   GraphTriangleCount[g]
 *   LocalClusteringCoefficient[g] / [g, v]
 *   GlobalClusteringCoefficient[g]
 *   MeanClusteringCoefficient[g]
 *
 * All four are exact (Integer / Rational), ignore EdgeWeight, and leave a
 * MIXED graph unevaluated -- Wolfram refuses mixed graphs for each of them.
 *
 * Undirected: t(v) = triangles through v, d(v) = degree.
 *   LocalClusteringCoefficient  t(v) / C(d(v), 2)   (0 when d(v) < 2)
 *   GlobalClusteringCoefficient 3 T / sum C(d(v), 2) (0 when that is 0)
 *   MeanClusteringCoefficient   mean of the local coefficients (0 for n = 0)
 *   GraphTriangleCount          T
 *
 * Directed (reverse-engineered from Mathematica 15, which documents none of
 * it; every case below was checked): a "triangle" is a directed 3-cycle
 * u->v->w->u. With c(v) the number of directed 3-cycles through v,
 *   local  c(v) / (in(v) out(v) - r(v)),  r(v) = neighbours joined to v in
 *          both directions (the number of in/out pairs that are the same
 *          vertex); 0 when that denominator is 0;
 *   global sum c(v) / sum (in(v) out(v) - r(v));
 *   mean   mean of the local coefficients;
 *   count  the number of directed 3-cycles (a transitive triple 1->2, 2->3,
 *          1->3 counts 0; a doubly-linked triangle counts 2).
 *
 * Algorithm: triangle listing over the underlying simple graph with the
 * degree-ordered orientation (each edge points from lower to higher (degree,
 * index) rank), which bounds every oriented out-degree by O(sqrt m) and gives
 * O(m^1.5) total work (Chiba-Nishizeki / Latapy). Each oriented edge carries a
 * two-bit record of which arcs exist, so the directed 3-cycle test is two
 * bit-ANDs per triangle. Rows are processed on the thread team with
 * per-thread mark arrays and counters.
 */

#include "graph_metrics.h"
#include "graph.h"
#include "expr.h"
#include "arithmetic.h"
#include "sym_names.h"
#include <gmp.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int n;
    int directed;
    int64_t* t;        /* triangles (undirected) / 3-cycles (directed) per vertex */
    int64_t* den;      /* local denominator per vertex                            */
    int64_t total;     /* T (undirected) or number of directed 3-cycles           */
} TriInfo;

static void tri_free(TriInfo* ti) {
    free(ti->t); free(ti->den);
    memset(ti, 0, sizeof(*ti));
}

/* Oriented adjacency: row v holds the neighbours of higher rank, each with its
 * arc flag from v's point of view (bit0: v->w exists, bit1: w->v exists). */
typedef struct {
    const int64_t* off;
    const int* adj;
    const unsigned char* fl;   /* NULL for an undirected graph (all flags 3)   */
    int n;
    int directed;
    int per_vertex;            /* also accumulate t[] (not needed for a count) */
    int64_t* t[32];
    int64_t tot[32];
    unsigned char* mark[32];
} TriCtx;

static void tri_work(void* vctx, int tid, int64_t lo, int64_t hi) {
    TriCtx* c = (TriCtx*)vctx;
    unsigned char* mark = c->mark[tid];
    int64_t* t = c->t[tid];
    int64_t tot = 0;
    const int64_t* off = c->off;
    const int* adj = c->adj;
    const unsigned char* fl = c->fl;
    for (int64_t u = lo; u < hi; u++) {
        int64_t a = off[u], b = off[u + 1];
        if (b - a < 2 && !c->per_vertex) continue;
        if (!c->directed) {
            for (int64_t p = a; p < b; p++) mark[adj[p]] = 1;
            int64_t tu = 0;
            for (int64_t p = a; p < b; p++) {
                int v = adj[p];
                int64_t tv = 0;
                for (int64_t q = off[v]; q < off[v + 1]; q++) {
                    int w = adj[q];
                    if (mark[w]) { tv++; if (c->per_vertex) t[w]++; }
                }
                tu += tv;
                if (c->per_vertex) t[v] += tv;
            }
            if (c->per_vertex) t[u] += tu;
            tot += tu;
        } else {
            for (int64_t p = a; p < b; p++) mark[adj[p]] = (unsigned char)(fl[p] + 1);
            for (int64_t p = a; p < b; p++) {
                int v = adj[p];
                int fuv = fl[p];
                for (int64_t q = off[v]; q < off[v + 1]; q++) {
                    int w = adj[q];
                    int m = mark[w];
                    if (!m) continue;
                    int fuw = m - 1, fvw = fl[q];
                    int cnt = ((fuv & 1) && (fvw & 1) && (fuw & 2))
                            + ((fuw & 1) && (fvw & 2) && (fuv & 2));
                    if (!cnt) continue;
                    t[u] += cnt; t[v] += cnt; t[w] += cnt;
                    tot += cnt;
                }
            }
        }
        for (int64_t p = a; p < b; p++) mark[adj[p]] = 0;
    }
    c->tot[tid] += tot;
}

typedef struct { int w; unsigned char f; } NbrFlag;

static int nbr_cmp(const void* a, const void* b) {
    int x = ((const NbrFlag*)a)->w, y = ((const NbrFlag*)b)->w;
    return (x > y) - (x < y);
}

/* Computes per-vertex triangle / 3-cycle counts and local denominators.
 * Returns 1 on success, 0 if g is invalid or mixed (or on allocation failure). */
static int tri_compute(const Expr* g, TriInfo* ti, int per_vertex) {
    memset(ti, 0, sizeof(*ti));
    int kind = gmet_graph_kind(g);
    if (kind < 0 || kind == GMET_KIND_MIXED) return 0;
    const int *eu, *ev;
    const unsigned char* edir;
    graph_edge_indices(g, &eu, &ev, &edir);
    int n = (int)g->data.function.args[0]->data.function.arg_count;
    size_t ne = g->data.function.args[1]->data.function.arg_count;
    size_t nn = (size_t)(n > 0 ? n : 1);
    int directed = kind == GMET_KIND_DIRECTED;
    ti->n = n; ti->directed = directed;
    ti->t = calloc(nn, sizeof(int64_t));
    ti->den = calloc(nn, sizeof(int64_t));
    if (!ti->t || !ti->den) { tri_free(ti); return 0; }
    if (n == 0) return 1;

    int64_t* ooff = NULL;
    int* oadj = NULL;
    unsigned char* ofl = NULL;
    int64_t no = 0;
    if (!directed) {
        /* Undirected: orient each edge straight from the endpoint arrays. */
        int* deg = calloc(nn, sizeof(int));
        int* rank = malloc(nn * sizeof(int));
        ooff = calloc(nn + 1, sizeof(int64_t));
        if (!deg || !rank || !ooff) { free(deg); free(rank); free(ooff); tri_free(ti); return 0; }
        /* One pass gives the degrees AND the out-degrees of the orientation
         * by index (each edge from its smaller to its larger endpoint). */
        for (size_t k = 0; k < ne; k++) {
            int a = eu[k], b = ev[k];
            deg[a]++; deg[b]++;
            ooff[(a < b ? a : b) + 1]++;
        }
        int maxd = 0;
        for (int v = 0; v < n; v++) {
            if (deg[v] > maxd) maxd = deg[v];
            ti->den[v] = (int64_t)deg[v] * (deg[v] - 1) / 2;
        }
        /* Orientation by index bounds the listing work by m * maxdeg; the
         * degree ordering bounds it by O(m^1.5) for any degree sequence. Use
         * the cheaper-to-build index order when maxdeg^2 <= 4m (bounded-degree
         * graphs), where the two bounds agree. */
        int by_index = (double)maxd * (double)maxd <= 4.0 * (double)ne;
        if (!by_index) {
            int* cnt = calloc((size_t)maxd + 2, sizeof(int));
            if (!cnt) { free(deg); free(rank); free(ooff); tri_free(ti); return 0; }
            for (int v = 0; v < n; v++) cnt[deg[v] + 1]++;
            for (int d = 0; d <= maxd; d++) cnt[d + 1] += cnt[d];
            for (int v = 0; v < n; v++) rank[v] = cnt[deg[v]]++;
            free(cnt);
            memset(ooff, 0, (nn + 1) * sizeof(int64_t));
            for (size_t k = 0; k < ne; k++) ooff[(rank[eu[k]] < rank[ev[k]] ? eu[k] : ev[k]) + 1]++;
        }
        free(deg);
        for (int v = 0; v < n; v++) ooff[v + 1] += ooff[v];
        no = ooff[n];
        oadj = malloc((size_t)(no + 1) * sizeof(int));
        int64_t* pos = malloc(nn * sizeof(int64_t));
        if (!oadj || !pos) { free(oadj); free(pos); free(rank); free(ooff); tri_free(ti); return 0; }
        memcpy(pos, ooff, nn * sizeof(int64_t));
        if (by_index) {
            for (size_t k = 0; k < ne; k++) {
                int a = eu[k], b = ev[k];
                if (a < b) oadj[pos[a]++] = b; else oadj[pos[b]++] = a;
            }
        } else {
            for (size_t k = 0; k < ne; k++) {
                int a = eu[k], b = ev[k];
                if (rank[a] < rank[b]) oadj[pos[a]++] = b; else oadj[pos[b]++] = a;
            }
        }
        free(pos); free(rank);
    } else {
    /* Directed: the underlying simple graph with arc flags, CSR by vertex. */
    int64_t* off = calloc(nn + 1, sizeof(int64_t));
    NbrFlag* nb = malloc((2 * ne + 1) * sizeof(NbrFlag));
    int64_t* fill = malloc((nn + 1) * sizeof(int64_t));
    int64_t* indeg = calloc(nn, sizeof(int64_t));
    int64_t* outdeg = calloc(nn, sizeof(int64_t));
    if (!off || !nb || !fill || !indeg || !outdeg) {
        free(off); free(nb); free(fill); free(indeg); free(outdeg); tri_free(ti); return 0;
    }
    for (size_t k = 0; k < ne; k++) { off[eu[k] + 1]++; off[ev[k] + 1]++; }
    for (int i = 0; i < n; i++) off[i + 1] += off[i];
    memcpy(fill, off, (nn + 1) * sizeof(int64_t));
    for (size_t k = 0; k < ne; k++) {
        int a = eu[k], b = ev[k];
        nb[fill[a]].w = b; nb[fill[a]].f = 1; fill[a]++;
        nb[fill[b]].w = a; nb[fill[b]].f = 2; fill[b]++;
        outdeg[a]++; indeg[b]++;
    }
    int* deg = malloc(nn * sizeof(int));
    if (!deg) { free(off); free(nb); free(fill); free(indeg); free(outdeg); tri_free(ti); return 0; }
    /* Merge u->v with v->u into one neighbour entry (sort + dedupe). */
    for (int v = 0; v < n; v++) {
        int64_t a = off[v], b = off[v + 1];
        if (b - a > 1) qsort(nb + a, (size_t)(b - a), sizeof(NbrFlag), nbr_cmp);
        int64_t w = a;
        int64_t recip = 0;
        for (int64_t p = a; p < b; p++) {
            if (w > a && nb[w - 1].w == nb[p].w) nb[w - 1].f |= nb[p].f;
            else nb[w++] = nb[p];
        }
        for (int64_t p = a; p < w; p++) if (nb[p].f == 3) recip++;
        fill[v] = w;                        /* end of the merged row */
        deg[v] = (int)(w - a);
        ti->den[v] = indeg[v] * outdeg[v] - recip;
    }
    free(indeg); free(outdeg);

    /* Rank by (degree, index) via counting sort on degree. */
    int maxd = 0;
    for (int v = 0; v < n; v++) if (deg[v] > maxd) maxd = deg[v];
    int* cnt = calloc((size_t)maxd + 2, sizeof(int));
    int* rank = malloc(nn * sizeof(int));
    ooff = calloc(nn + 1, sizeof(int64_t));
    if (!cnt || !rank || !ooff) {
        free(cnt); free(rank); free(ooff); free(deg); free(off); free(nb); free(fill);
        tri_free(ti); return 0;
    }
    for (int v = 0; v < n; v++) cnt[deg[v] + 1]++;
    for (int d = 0; d <= maxd; d++) cnt[d + 1] += cnt[d];
    for (int v = 0; v < n; v++) rank[v] = cnt[deg[v]]++;
    free(cnt);
    for (int v = 0; v < n; v++) {
        for (int64_t p = off[v]; p < fill[v]; p++) if (rank[nb[p].w] > rank[v]) ooff[v + 1]++;
        no += ooff[v + 1];
    }
    for (int v = 0; v < n; v++) ooff[v + 1] += ooff[v];
    oadj = malloc((size_t)(no + 1) * sizeof(int));
    ofl = malloc((size_t)(no + 1));
    if (!oadj || !ofl) {
        free(oadj); free(ofl); free(rank); free(ooff); free(deg); free(off); free(nb); free(fill);
        tri_free(ti); return 0;
    }
    for (int v = 0; v < n; v++) {
        int64_t q = ooff[v];
        for (int64_t p = off[v]; p < fill[v]; p++)
            if (rank[nb[p].w] > rank[v]) { oadj[q] = nb[p].w; ofl[q] = nb[p].f; q++; }
    }
    free(rank); free(deg); free(off); free(nb); free(fill);
    }

    TriCtx c;
    memset(&c, 0, sizeof(c));
    c.off = ooff; c.adj = oadj; c.fl = ofl; c.n = n; c.directed = directed;
    c.per_vertex = per_vertex || directed;
    int nt = gmet_cap_threads(gmet_thread_count((double)no * 16.0), 9.0 * (double)n);
    int ok = 1;
    for (int t = 0; t < nt; t++) {
        c.mark[t] = calloc(nn, 1);
        c.t[t] = t == 0 ? ti->t : calloc(nn, sizeof(int64_t));
        if (!c.mark[t] || !c.t[t]) ok = 0;
    }
    if (ok) {
        gmet_parallel_for(n, 256, nt, tri_work, &c);
        for (int t = 1; t < nt; t++)
            for (int v = 0; v < n; v++) ti->t[v] += c.t[t][v];
        for (int t = 0; t < nt; t++) ti->total += c.tot[t];
    }
    for (int t = 0; t < nt; t++) { free(c.mark[t]); if (t > 0) free(c.t[t]); }
    free(ooff); free(oadj); free(ofl);
    if (!ok) { tri_free(ti); return 0; }
    return 1;
}

/* ---- Builtins -------------------------------------------------------------- */

Expr* builtin_graph_triangle_count(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    TriInfo ti;
    if (!tri_compute(res->data.function.args[0], &ti, 0)) return NULL;
    Expr* out = expr_new_integer(ti.total);
    tri_free(&ti);
    return out;
}

static Expr* local_value(const TriInfo* ti, int v) {
    if (ti->den[v] <= 0 || ti->t[v] == 0) return expr_new_integer(0);
    return make_rational(ti->t[v], ti->den[v]);
}

Expr* builtin_local_clustering_coefficient(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    int vi = -1;
    if (argc == 2) {
        vi = graph_vertex_position(g, res->data.function.args[1]);
        if (vi < 0) return NULL;
    }
    Expr* cached = argc == 1 ? gmet_cache_get("LocalClusteringCoefficient", res) : NULL;
    if (cached) return cached;
    TriInfo ti;
    if (!tri_compute(g, &ti, 1)) return NULL;
    Expr* out;
    if (argc == 2) out = local_value(&ti, vi);
    else {
        /* Degree sequences repeat, so most coefficients repeat: build each
         * distinct (t, den) value once and share it by reference. */
        Expr** it = malloc((size_t)(ti.n > 0 ? ti.n : 1) * sizeof(Expr*));
        size_t cap = 64;
        while (cap < (size_t)ti.n * 2) cap <<= 1;
        int* slot = malloc(cap * sizeof(int));
        if (slot) for (size_t i = 0; i < cap; i++) slot[i] = -1;
        for (int v = 0; v < ti.n; v++) {
            if (!slot) { it[v] = local_value(&ti, v); continue; }
            uint64_t key = (uint64_t)ti.t[v] * 0x9E3779B97F4A7C15ULL ^ (uint64_t)ti.den[v] * 0xC2B2AE3D27D4EB4FULL;
            size_t h = (size_t)(key >> 17) & (cap - 1);
            while (slot[h] >= 0 && !(ti.t[slot[h]] == ti.t[v] && ti.den[slot[h]] == ti.den[v]))
                h = (h + 1) & (cap - 1);
            if (slot[h] >= 0) it[v] = expr_copy(it[slot[h]]);
            else { slot[h] = v; it[v] = local_value(&ti, v); }
        }
        free(slot);
        out = expr_new_function(expr_new_symbol(SYM_List), it, (size_t)ti.n);
        free(it);
        gmet_cache_put("LocalClusteringCoefficient", res, out);
    }
    tri_free(&ti);
    return out;
}

Expr* builtin_global_clustering_coefficient(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    TriInfo ti;
    if (!tri_compute(res->data.function.args[0], &ti, 1)) return NULL;
    int64_t num = 0, den = 0;
    for (int v = 0; v < ti.n; v++) { num += ti.t[v]; den += ti.den[v]; }
    tri_free(&ti);
    if (den <= 0 || num == 0) return expr_new_integer(0);
    return make_rational(num, den);
}

typedef struct { int64_t den, num; } DenNum;

static int dn_cmp(const void* a, const void* b) {
    int64_t x = ((const DenNum*)a)->den, y = ((const DenNum*)b)->den;
    return (x > y) - (x < y);
}

Expr* builtin_mean_clustering_coefficient(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    TriInfo ti;
    if (!tri_compute(res->data.function.args[0], &ti, 1)) return NULL;
    int n = ti.n;
    if (n == 0) { tri_free(&ti); return expr_new_integer(0); }
    /* Exact mean: group numerators by denominator (few distinct values), then
     * add the groups as GMP rationals. */
    DenNum* dn = malloc((size_t)n * sizeof(DenNum));
    if (!dn) { tri_free(&ti); return NULL; }
    int k = 0;
    for (int v = 0; v < n; v++)
        if (ti.den[v] > 0 && ti.t[v] > 0) { dn[k].den = ti.den[v]; dn[k].num = ti.t[v]; k++; }
    tri_free(&ti);
    qsort(dn, (size_t)k, sizeof(DenNum), dn_cmp);
    mpq_t acc, term;
    mpq_init(acc); mpq_init(term);
    mpz_t z;
    mpz_init(z);
    for (int i = 0; i < k;) {
        int64_t d = dn[i].den;
        int64_t s = 0;
        while (i < k && dn[i].den == d) s += dn[i++].num;
        mpz_set_si(z, (long)s);
        mpq_set_z(term, z);
        mpz_set_si(z, (long)d);
        mpq_set_den(term, z);
        mpq_canonicalize(term);
        mpq_add(acc, acc, term);
    }
    free(dn);
    mpz_set_si(z, (long)n);
    mpq_set_z(term, z);
    mpq_div(acc, acc, term);
    Expr* out = make_rational_mpz(mpq_numref(acc), mpq_denref(acc));
    mpq_clear(acc); mpq_clear(term); mpz_clear(z);
    return out;
}

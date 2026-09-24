/* galg_flow.c - network flows and cuts.
 *
 *   FindMaximumFlow[g, s, t]            the maximum flow value from s to t
 *   FindMaximumFlow[g, s, t, "prop"]    "FlowValue", "FlowMatrix", "EdgeList"
 *   FindMinimumCut[g]                   {value, {part1, part2}}: a global
 *                                       minimum edge cut
 *   FindEdgeCut[g] / FindEdgeCut[g,s,t] the edges of a minimum (s-t) edge cut
 *   FindVertexCut[g] / [g, s, t]        a minimum (s-t) vertex separator
 *   EdgeConnectivity[g] / [g, s, t]     the (s-t) edge connectivity
 *
 * Capacities (Wolfram-faithful, verified against Mathematica 15):
 *   - FindMaximumFlow takes capacities from its EdgeCapacity -> {c1, ...}
 *     option (EdgeList order) and IGNORES EdgeWeight; without the option every
 *     edge has capacity 1. VertexCapacity -> {c1, ...} (VertexList order) caps
 *     the flow through each vertex other than the sources and sinks.
 *   - The cut family (FindMinimumCut, FindEdgeCut, EdgeConnectivity) uses
 *     EdgeWeight as the capacity when g carries one, else 1.
 *   - Integer capacities give exact Integer results. Rational or Real ones give
 *     a Real, as Mathematica does. Infinity is an allowed capacity (a flow of
 *     Infinity is reported as Infinity -- Mathematica's answer there is not
 *     meaningful). A negative or symbolic capacity leaves the call unevaluated.
 *   - An UndirectedEdge carries flow either way up to its capacity; a
 *     DirectedEdge only forwards. Mixed graphs are fine.
 *   - FindVertexCut works on the underlying undirected graph (as Mathematica
 *     does): FindVertexCut[g] of a complete graph is its first n-1 vertices
 *     (but {} when g has a directed edge and n >= 3, again as Mathematica),
 *     and adjacent s, t give {}.
 *
 * Exactness: all arithmetic is on int64. Non-integer capacities are scaled by
 * a common power of two -- chosen so every sum fits -- which makes machine
 * reals exact dyadic integers (their exact max flow is then computed), and a
 * Rational goes through its double value, which is precisely the precision of
 * the Real Mathematica returns for it.
 *
 * Algorithms:
 *   - Max flow: Dinic's algorithm on a CSR residual network (arcs of one edge
 *     are adjacent reverse pairs), BFS levels truncated at the sink's level,
 *     and an iterative blocking-flow search with current-arc pointers. On unit
 *     capacities this is O(E sqrt(E)); in practice a handful of phases.
 *   - Global min cut of an undirected graph: Nagamochi-Ibaraki. Each round
 *     computes a maximum-adjacency order, lowers the bound lambda to the least
 *     weighted degree, and contracts EVERY edge whose MA-order attachment
 *     q(e) >= lambda (those edges carry lambda(u,v) >= lambda, so contracting
 *     them cannot destroy a lighter cut) plus the last two vertices of the
 *     order. It returns exactly Stoer-Wagner's answer but typically contracts
 *     most of the graph per round instead of one pair.
 *   - Global min cut of a (partly) directed graph: min over v != v0 of the
 *     flows v0 -> v and v -> v0, each capped at the running bound.
 *   - Vertex connectivity / separators: Even's split-vertex construction
 *     (v_in -> v_out with capacity 1) and the Esfahanian-Hakimi pair
 *     selection: a minimum-degree vertex against each non-neighbour, then
 *     each non-adjacent pair of its neighbours -- every flow capped at the
 *     running bound.
 *
 * Tie-breaking for cut sides (the answer is unique only up to ties): an s-t
 * EDGE cut is the one closest to s, an s-t VERTEX cut the one closest to t --
 * Mathematica's choices on every probe. Cut edges come in EdgeList order,
 * vertex sets in VertexList order. FindMinimumCut lists the source side first
 * for a directed graph, and the side without VertexList[g][[1]] first for an
 * undirected one.
 *
 * Memory (SPEC section 4): fresh results; scratch freed on every path.
 */

#include "graph.h"
#include "graph_algos.h"
#include "expr.h"
#include "sym_names.h"
#include "pack.h"
#include "ndarray.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ======================================================================== *
 * Capacities
 * ======================================================================== */

typedef struct {
    int64_t* c;     /* scaled capacity per edge; GF_INF_MARK for Infinity      */
    int      real;  /* 1 if the results are Real                               */
    int      shift; /* true value = c / 2^shift                                */
    int      any_inf;
    int64_t  inf;   /* the finite stand-in for Infinity (sum of finite caps+1) */
} GfCap;

#define GF_INF_MARK INT64_MIN
/* Every scaled sum must stay below this: one flow value, or the finite
 * stand-in for Infinity, is at most (sum of caps + 1). */
#define GF_SUM_LIMIT ((int64_t)1 << 61)

static void gf_cap_free(GfCap* c) { free(c->c); c->c = NULL; }

/* Parse a capacity list of length m (NULL: all 1). Returns 1, or 0 when some
 * entry is not a usable (non-negative numeric or Infinity) capacity or the
 * list is malformed -- the caller stays unevaluated. */
static int gf_caps_parse(const Expr* lst_in, long m, GfCap* out) {
    memset(out, 0, sizeof(*out));
    out->c = malloc((size_t)(m > 0 ? m : 1) * sizeof(int64_t));
    if (!out->c) return 0;
    if (!lst_in) {
        for (long k = 0; k < m; k++) out->c[k] = 1;
        out->inf = m + 1;
        return 1;
    }
    double* d = malloc((size_t)(m > 0 ? m : 1) * sizeof(double));
    if (!d) { gf_cap_free(out); return 0; }
    int ok = 1, all_int = 1;
    /* Packed fast path: read the machine buffer, no unpacking. */
    if (lst_in->type == EXPR_NDARRAY && is_packed_list(lst_in)
        && lst_in->data.ndarray.rank == 1 && lst_in->data.ndarray.dims[0] == m
        && (lst_in->data.ndarray.dtype == NDT_INT64 || lst_in->data.ndarray.dtype == NDT_FLOAT64)) {
        if (lst_in->data.ndarray.dtype == NDT_INT64) {
            const int64_t* src = (const int64_t*)lst_in->data.ndarray.data;
            for (long k = 0; k < m && ok; k++) {
                if (src[k] < 0) ok = 0;
                out->c[k] = src[k]; d[k] = (double)src[k];
            }
        } else {
            const double* src = (const double*)lst_in->data.ndarray.data;
            all_int = 0;
            for (long k = 0; k < m && ok; k++) {
                if (!(src[k] >= 0) || isinf(src[k])) ok = 0;
                d[k] = src[k]; out->c[k] = 0;
            }
        }
        goto parsed;
    }
    Expr* owned = NULL;
    const Expr* lst = galg_plain_list(lst_in, &owned);
    if (!lst || (long)lst->data.function.arg_count != m) { expr_free(owned); free(d); gf_cap_free(out); return 0; }
    for (long k = 0; k < m && ok; k++) {
        const Expr* w = lst->data.function.args[k];
        if (w->type == EXPR_INTEGER) {
            if (w->data.integer < 0) ok = 0;
            d[k] = (double)w->data.integer;
            out->c[k] = w->data.integer;
        } else if (galg_is_symbol(w, "Infinity")) {
            out->c[k] = GF_INF_MARK; d[k] = -1; out->any_inf = 1;
        } else {
            double x = graph_weight_to_double(w);
            if (!(x >= 0) || isinf(x)) ok = 0;     /* NaN, negative, overflow */
            d[k] = x; all_int = 0;
        }
    }
    expr_free(owned);
parsed:
    if (!ok) { free(d); gf_cap_free(out); return 0; }

    if (all_int) {
        /* Exact integers: only the sum bound to check. */
        int64_t sum = 0;
        for (long k = 0; k < m; k++) {
            if (out->c[k] == GF_INF_MARK) continue;
            if (out->c[k] > GF_SUM_LIMIT - sum) { ok = 0; break; }
            sum += out->c[k];
        }
        free(d);
        if (!ok) { gf_cap_free(out); return 0; }
        out->inf = sum + 1;
    } else {
        /* Reals: scale by 2^shift, the smallest power making every capacity
         * an integer, lowered as needed so the total stays under the limit
         * (then capacities are rounded -- far below machine precision of the
         * Real result). */
        out->real = 1;
        double maxd = 0, sumd = 0;
        int need = 0;
        for (long k = 0; k < m; k++) {
            if (out->c[k] == GF_INF_MARK || d[k] == 0) continue;
            if (d[k] > maxd) maxd = d[k];
            sumd += d[k];
            int e;
            double f = frexp(d[k], &e);          /* d = f * 2^e, f in [0.5,1) */
            /* mantissa bits actually used */
            double mant = ldexp(f, 53);
            int tz = 0;
            while (tz < 53 && fmod(mant, 2.0) == 0) { mant /= 2; tz++; }
            int frac = 53 - tz - e;               /* bits below the binary point */
            if (frac > need) need = frac;
        }
        int shift = need;
        if (sumd > 0) {
            /* largest shift with sumd * 2^shift < GF_SUM_LIMIT / 2 */
            int e;
            frexp(sumd, &e);                      /* sumd < 2^e */
            int maxshift = 59 - e;
            if (shift > maxshift) shift = maxshift;
        }
        out->shift = shift;
        int64_t sum = 0;
        for (long k = 0; k < m; k++) {
            if (out->c[k] == GF_INF_MARK) continue;
            out->c[k] = (int64_t)llround(ldexp(d[k], shift));
            sum += out->c[k];
        }
        free(d);
        if (sum >= GF_SUM_LIMIT) { gf_cap_free(out); return 0; }
        out->inf = sum + 1;
    }
    for (long k = 0; k < m; k++) if (out->c[k] == GF_INF_MARK) out->c[k] = out->inf;
    return 1;
}

/* Parse EdgeCapacity (el, m entries; NULL: all 1) and VertexCapacity (vl, n
 * entries) TOGETHER, so both lists share one exactness flag and one scale
 * 2^shift -- they meet in a single network, where capacities on different
 * scales would be summed wrongly. The two halves of one joint parse are split
 * into ec / vc, each keeping its own any_inf; both carry the joint inf, a
 * finite stand-in exceeding every finite sum of either. Returns 1, or 0 as
 * gf_caps_parse does. */
static int gf_caps_parse_joint(const Expr* el, long m, const Expr* vl, long n,
                               GfCap* ec, GfCap* vc) {
    memset(ec, 0, sizeof(*ec));
    memset(vc, 0, sizeof(*vc));
    Expr *eo = NULL, *vo = NULL;
    const Expr* ep = el ? galg_plain_list(el, &eo) : NULL;
    const Expr* vp = galg_plain_list(vl, &vo);
    int ok = (!el || (ep && (long)ep->data.function.arg_count == m))
             && vp && (long)vp->data.function.arg_count == n;
    GfCap all;
    memset(&all, 0, sizeof(all));
    if (ok) {
        Expr** items = malloc((size_t)(m + n > 0 ? m + n : 1) * sizeof(Expr*));
        ok = items != NULL;
        if (ok) {
            for (long k = 0; k < m; k++)
                items[k] = ep ? expr_copy(ep->data.function.args[k]) : expr_new_integer(1);
            for (long k = 0; k < n; k++) items[m + k] = expr_copy(vp->data.function.args[k]);
            Expr* both = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)(m + n));
            free(items);
            ok = gf_caps_parse(both, m + n, &all);
            expr_free(both);
        }
    }
    expr_free(eo); expr_free(vo);
    if (!ok) return 0;
    ec->c = malloc((size_t)(m > 0 ? m : 1) * sizeof(int64_t));
    vc->c = malloc((size_t)(n > 0 ? n : 1) * sizeof(int64_t));
    if (!ec->c || !vc->c) { gf_cap_free(ec); gf_cap_free(vc); gf_cap_free(&all); return 0; }
    memcpy(ec->c, all.c, (size_t)m * sizeof(int64_t));
    memcpy(vc->c, all.c + m, (size_t)n * sizeof(int64_t));
    ec->real = vc->real = all.real;
    ec->shift = vc->shift = all.shift;
    ec->inf = vc->inf = all.inf;
    /* an Infinity entry was replaced by all.inf, which no finite entry reaches */
    for (long k = 0; k < m; k++) if (ec->c[k] == all.inf) ec->any_inf = 1;
    for (long k = 0; k < n; k++) if (vc->c[k] == all.inf) vc->any_inf = 1;
    gf_cap_free(&all);
    return 1;
}

/* The Expr for a scaled flow/cut value. */
static Expr* gf_value(const GfCap* c, int64_t v) {
    if (c->any_inf && v >= c->inf) return expr_new_symbol(SYM_Infinity);
    if (!c->real) return expr_new_integer(v);
    return expr_new_real(ldexp((double)v, -c->shift));
}

/* ======================================================================== *
 * Residual network and Dinic
 * ======================================================================== */

typedef struct {
    int      n;
    long     na;
    long*    first;   /* n + 1: arcs of v are first[v] .. first[v+1]-1         */
    int*     to;
    long*    rev;     /* index of the reverse arc                              */
    int64_t* res;     /* residual capacity                                     */
    int64_t* cap0;    /* initial residual (for resets), may be NULL            */
    /* Dinic scratch */
    int*     level;
    long*    it;
    int*     queue;
    long*    pstack;
} GfNet;

static void gf_net_free(GfNet* N) {
    if (!N) return;
    free(N->first); free(N->to); free(N->rev); free(N->res); free(N->cap0);
    free(N->level); free(N->it); free(N->queue); free(N->pstack);
    free(N);
}

/* Build from np arc PAIRS: pair p is the arc pu[p] -> pv[p] with capacity
 * cf[p] and its reverse with capacity cb[p] (cb = cf for an undirected edge,
 * 0 for a directed one). pos_out[p] (when non-NULL) receives the CSR index of
 * the forward arc. keep_cap0 stores the initial residuals for gf_net_reset. */
static GfNet* gf_net_build(int n, long np, const int* pu, const int* pv,
                           const int64_t* cf, const int64_t* cb, long* pos_out,
                           int keep_cap0) {
    GfNet* N = calloc(1, sizeof(GfNet));
    if (!N) return NULL;
    size_t nn = (size_t)n + 1, na = (size_t)(2 * np) + 1;
    N->n = n; N->na = 2 * np;
    N->first = calloc(nn + 1, sizeof(long));
    N->to    = malloc(na * sizeof(int));
    N->rev   = malloc(na * sizeof(long));
    N->res   = malloc(na * sizeof(int64_t));
    N->level = malloc(nn * sizeof(int));
    N->it    = malloc(nn * sizeof(long));
    N->queue = malloc(nn * sizeof(int));
    N->pstack = malloc(nn * sizeof(long));
    if (keep_cap0) N->cap0 = malloc(na * sizeof(int64_t));
    if (!N->first || !N->to || !N->rev || !N->res || !N->level || !N->it
        || !N->queue || !N->pstack || (keep_cap0 && !N->cap0)) { gf_net_free(N); return NULL; }
    for (long p = 0; p < np; p++) { N->first[pu[p] + 1]++; N->first[pv[p] + 1]++; }
    for (int v = 0; v < n; v++) N->first[v + 1] += N->first[v];
    long* fill = malloc(nn * sizeof(long));
    if (!fill) { gf_net_free(N); return NULL; }
    for (int v = 0; v < n; v++) fill[v] = N->first[v];
    for (long p = 0; p < np; p++) {
        long a = fill[pu[p]]++, b = fill[pv[p]]++;
        N->to[a] = pv[p]; N->to[b] = pu[p];
        N->rev[a] = b; N->rev[b] = a;
        N->res[a] = cf[p]; N->res[b] = cb[p];
        if (pos_out) pos_out[p] = a;
    }
    free(fill);
    if (keep_cap0) memcpy(N->cap0, N->res, (size_t)N->na * sizeof(int64_t));
    return N;
}

static void gf_net_reset(GfNet* N) {
    memcpy(N->res, N->cap0, (size_t)N->na * sizeof(int64_t));
}

/* BFS levels from s over positive residual arcs; stops once t's level is
 * settled. Returns 1 iff t is reachable. */
static int gf_bfs(GfNet* N, int s, int t) {
    int n = N->n;
    for (int v = 0; v < n; v++) N->level[v] = -1;
    int head = 0, tail = 0;
    N->level[s] = 0; N->queue[tail++] = s;
    while (head < tail) {
        int v = N->queue[head++];
        if (N->level[t] >= 0 && N->level[v] >= N->level[t]) break;
        int lv = N->level[v] + 1;
        for (long e = N->first[v]; e < N->first[v + 1]; e++) {
            int w = N->to[e];
            if (N->res[e] > 0 && N->level[w] < 0) {
                N->level[w] = lv;
                N->queue[tail++] = w;
            }
        }
    }
    return N->level[t] >= 0;
}

/* Dinic from s to t, stopping once the flow reaches `limit`. Returns the flow
 * pushed by this call (the residuals keep it). */
static int64_t gf_dinic(GfNet* N, int s, int t, int64_t limit) {
    int64_t flow = 0;
    if (s == t) return 0;
    long polls = 0;
    while (flow < limit && gf_bfs(N, s, t)) {
        if ((++polls & 63) == 0) galg_poll();
        for (int v = 0; v < N->n; v++) N->it[v] = N->first[v];
        int top = 0, v = s;
        for (;;) {
            if (v == t) {
                int64_t f = limit - flow;
                for (int i = 0; i < top; i++) if (N->res[N->pstack[i]] < f) f = N->res[N->pstack[i]];
                int cut = -1;
                for (int i = 0; i < top; i++) {
                    long e = N->pstack[i];
                    N->res[e] -= f; N->res[N->rev[e]] += f;
                    if (cut < 0 && N->res[e] == 0) cut = i;
                }
                flow += f;
                if (flow >= limit) return flow;
                top = cut;                           /* retreat to the bottleneck */
                v = top == 0 ? s : N->to[N->pstack[top - 1]];
                continue;
            }
            long e = N->it[v], end = N->first[v + 1];
            int lv = N->level[v] + 1;
            while (e < end && !(N->res[e] > 0 && N->level[N->to[e]] == lv)) e++;
            N->it[v] = e;
            if (e == end) {                          /* dead end: prune v */
                N->level[v] = -1;
                if (top == 0) break;
                top--;
                v = top == 0 ? s : N->to[N->pstack[top - 1]];
                N->it[v]++;
                continue;
            }
            N->pstack[top++] = e;
            v = N->to[e];
        }
    }
    return flow;
}

/* mark[v] = 1 for v reachable from s in the residual network. */
static void gf_reach_from(const GfNet* N, int s, char* mark) {
    memset(mark, 0, (size_t)N->n);
    int head = 0, tail = 0;
    mark[s] = 1; N->queue[tail++] = s;
    while (head < tail) {
        int v = N->queue[head++];
        for (long e = N->first[v]; e < N->first[v + 1]; e++) {
            int w = N->to[e];
            if (N->res[e] > 0 && !mark[w]) { mark[w] = 1; N->queue[tail++] = w; }
        }
    }
}

/* mark[v] = 1 for v that can reach t in the residual network. */
static void gf_reach_to(const GfNet* N, int t, char* mark) {
    memset(mark, 0, (size_t)N->n);
    int head = 0, tail = 0;
    mark[t] = 1; N->queue[tail++] = t;
    while (head < tail) {
        int w = N->queue[head++];
        for (long e = N->first[w]; e < N->first[w + 1]; e++) {
            int v = N->to[e];                       /* arc w->v; v->w is rev */
            if (N->res[N->rev[e]] > 0 && !mark[v]) { mark[v] = 1; N->queue[tail++] = v; }
        }
    }
}

/* ======================================================================== *
 * Edge networks of a graph
 * ======================================================================== */

/* One arc pair per edge of g: forward capacity cap[k], reverse cap[k] when
 * undirected, 0 when directed. pos[k] = CSR index of edge k's forward arc.
 * Extra pairs (super source/sink hookups) may be appended by the caller via
 * xu/xv/xc (nx of them, all directed). */
static GfNet* gf_edge_net(const Expr* g, const int64_t* cap, long* pos, int nextra_v,
                          long nx, const int* xu, const int* xv, const int64_t* xc,
                          int keep_cap0) {
    const int *eu, *ev;
    const unsigned char* edir;
    if (!graph_edge_indices(g, &eu, &ev, &edir)) return NULL;
    int n = galg_nv(g);
    long m = galg_ne(g), np = m + nx;
    int* pu = malloc((size_t)(np > 0 ? np : 1) * sizeof(int));
    int* pv = malloc((size_t)(np > 0 ? np : 1) * sizeof(int));
    int64_t* cf = malloc((size_t)(np > 0 ? np : 1) * sizeof(int64_t));
    int64_t* cb = malloc((size_t)(np > 0 ? np : 1) * sizeof(int64_t));
    long* ppos = malloc((size_t)(np > 0 ? np : 1) * sizeof(long));
    GfNet* N = NULL;
    if (pu && pv && cf && cb && ppos) {
        for (long k = 0; k < m; k++) {
            pu[k] = eu[k]; pv[k] = ev[k];
            cf[k] = cap[k]; cb[k] = edir[k] ? 0 : cap[k];
        }
        for (long k = 0; k < nx; k++) {
            pu[m + k] = xu[k]; pv[m + k] = xv[k]; cf[m + k] = xc[k]; cb[m + k] = 0;
        }
        N = gf_net_build(n + nextra_v, np, pu, pv, cf, cb, ppos, keep_cap0);
        if (N && pos) memcpy(pos, ppos, (size_t)(m > 0 ? m : 1) * sizeof(long));
    }
    free(pu); free(pv); free(cf); free(cb); free(ppos);
    return N;
}

/* Resolve a vertex-or-list-of-vertices argument into indices. Returns the
 * count (>= 1), 0 if any element is not a vertex. *idx is malloc'd. */
static int gf_terminals(const Expr* g, const Expr* a, int** idx) {
    *idx = NULL;
    int p = galg_vertex_arg(g, a);
    if (p >= 0) {
        *idx = malloc(sizeof(int));
        if (!*idx) return 0;
        (*idx)[0] = p;
        return 1;
    }
    if (!graph_is_list(a) || a->data.function.arg_count == 0) return 0;
    int k = (int)a->data.function.arg_count;
    *idx = malloc((size_t)k * sizeof(int));
    if (!*idx) return 0;
    for (int i = 0; i < k; i++) {
        (*idx)[i] = galg_vertex_arg(g, a->data.function.args[i]);
        if ((*idx)[i] < 0) { free(*idx); *idx = NULL; return 0; }
    }
    return k;
}

/* ======================================================================== *
 * FindMaximumFlow
 * ======================================================================== */

enum { GF_PROP_VALUE, GF_PROP_MATRIX, GF_PROP_EDGES };

/* Dense n x n matrix of flows (packed when uniform). */
static Expr* gf_flow_matrix(int n, const int* tail, const int* head, const int64_t* f,
                            long cnt, const GfCap* c) {
    size_t nn = (size_t)n;
    if (nn * nn > (size_t)64 * 1024 * 1024) return NULL;   /* refuse > 64M cells */
    Expr** rows = malloc((nn > 0 ? nn : 1) * sizeof(Expr*));
    if (!rows) return NULL;
    int64_t* acc = calloc((nn > 0 ? nn * nn : 1), sizeof(int64_t));
    if (!acc) { free(rows); return NULL; }
    for (long k = 0; k < cnt; k++) acc[(size_t)tail[k] * nn + (size_t)head[k]] += f[k];
    Expr** cell = malloc((nn > 0 ? nn : 1) * sizeof(Expr*));
    if (!cell) { free(rows); free(acc); return NULL; }
    for (size_t i = 0; i < nn; i++) {
        for (size_t j = 0; j < nn; j++) cell[j] = gf_value(c, acc[i * nn + j]);
        rows[i] = pack_offer(expr_new_function(expr_new_symbol(SYM_List), cell, nn));
    }
    Expr* out = pack_offer(expr_new_function(expr_new_symbol(SYM_List), rows, nn));
    free(cell); free(rows); free(acc);
    return out;
}

Expr* builtin_find_maximum_flow(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 3) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    size_t nopt = galg_first_option(res, 3);
    int prop = GF_PROP_VALUE;
    if (nopt == 4) {
        const Expr* p = res->data.function.args[3];
        if (galg_is_string(p, "FlowValue")) prop = GF_PROP_VALUE;
        else if (galg_is_string(p, "FlowMatrix")) prop = GF_PROP_MATRIX;
        else if (galg_is_string(p, "EdgeList")) prop = GF_PROP_EDGES;
        else return NULL;
    } else if (nopt != 3) return NULL;
    static const char* const names[] = { "EdgeCapacity", "VertexCapacity", NULL };
    const Expr* ov[2];
    if (!galg_split_options(res, nopt, names, ov)) return NULL;
    if (ov[0] && galg_is_symbol(ov[0], "Automatic")) ov[0] = NULL;
    if (ov[1] && galg_is_symbol(ov[1], "Automatic")) ov[1] = NULL;

    int n = galg_nv(g);
    long m = galg_ne(g);
    int *S = NULL, *T = NULL;
    int ns = gf_terminals(g, res->data.function.args[1], &S);
    int nt = ns ? gf_terminals(g, res->data.function.args[2], &T) : 0;
    if (!ns || !nt) { free(S); free(T); return NULL; }

    GfCap cap, vcap;
    memset(&vcap, 0, sizeof(vcap));
    int capok = ov[1] ? gf_caps_parse_joint(ov[0], m, ov[1], n, &cap, &vcap)
                      : gf_caps_parse(ov[0], m, &cap);
    if (!capok) { free(S); free(T); return NULL; }

    /* A terminal in both S and T: Mathematica answers 0 for s == t. */
    char* role = calloc((size_t)n + 1, 1);     /* 1 source, 2 sink */
    if (!role) { gf_cap_free(&cap); gf_cap_free(&vcap); free(S); free(T); return NULL; }
    int overlap = 0;
    for (int i = 0; i < ns; i++) role[S[i]] |= 1;
    for (int i = 0; i < nt; i++) { if (role[T[i]] & 1) overlap = 1; role[T[i]] |= 2; }

    Expr* out = NULL;
    const int *eu, *ev;
    const unsigned char* edir;
    graph_edge_indices(g, &eu, &ev, &edir);

    /* Network layout: vertices 0..n-1 (the "in" copy when vertex capacities
     * are used; the "out" copy is n + v), then the super source/sink when
     * there are several terminals. */
    int split = ov[1] != NULL;
    int base = split ? 2 * n : n;
    int multi = (ns > 1 || nt > 1);
    int src = multi ? base : S[0], snk = multi ? base + 1 : T[0];
    /* under splitting, flow enters s at its in copy and leaves t from its out
     * copy, so the terminals' own capacity arcs are on every path */
    if (split && !multi) snk = n + T[0];
    int nvx = base + (multi ? 2 : 0);
    long mu = 0;                                 /* undirected edges */
    for (long k = 0; k < m; k++) if (!edir[k]) mu++;
    /* pairs: one per edge; under splitting an undirected edge needs a second
     * directed pair (b_out -> a_in), plus the n internal v_in -> v_out arcs;
     * plus the super-terminal hookups */
    long np = m + (split ? mu + n : 0) + (multi ? ns + nt : 0);
    int* pu = malloc((size_t)(np > 0 ? np : 1) * sizeof(int));
    int* pv = malloc((size_t)(np > 0 ? np : 1) * sizeof(int));
    int64_t* cf = malloc((size_t)(np > 0 ? np : 1) * sizeof(int64_t));
    int64_t* cb = malloc((size_t)(np > 0 ? np : 1) * sizeof(int64_t));
    long* pos = malloc((size_t)(np > 0 ? np : 1) * sizeof(long));
    GfNet* N = NULL;
    /* With VertexCapacity both lists come from one joint parse, so cap.inf ==
     * vcap.inf is the single finite stand-in for Infinity: every finite cut is
     * below it, every cut through Infinity capacities reaches it. */
    int64_t total = cap.inf;
    if (!pu || !pv || !cf || !cb || !pos) goto done;
    {
        long p = 0;
        for (long k = 0; k < m; k++, p++) {
            /* with vertex splitting, an edge a->b runs a_out -> b_in; without
             * it one pair with symmetric capacity models an undirected edge */
            pu[p] = split ? n + eu[k] : eu[k]; pv[p] = ev[k];
            cf[p] = cap.c[k];
            cb[p] = (edir[k] || split) ? 0 : cap.c[k];
        }
        if (split) {
            for (long k = 0; k < m; k++) if (!edir[k]) {
                pu[p] = n + ev[k]; pv[p] = eu[k]; cf[p] = cap.c[k]; cb[p] = 0; p++;
            }
            for (int v = 0; v < n; v++, p++) {
                pu[p] = v; pv[p] = n + v;
                /* terminals are capped too, as in Mathematica
                 * (VertexCapacity -> {2, 10, 10} on 1-2-3 gives 2) */
                cf[p] = vcap.c[v];
                cb[p] = 0;
            }
        }
        if (multi) {
            for (int i = 0; i < ns; i++, p++) {
                pu[p] = src; pv[p] = S[i]; cf[p] = total; cb[p] = 0;
            }
            for (int i = 0; i < nt; i++, p++) {
                pu[p] = split ? n + T[i] : T[i]; pv[p] = snk; cf[p] = total; cb[p] = 0;
            }
        }
    }
    N = gf_net_build(nvx, np, pu, pv, cf, cb, pos, 0);
    if (!N) goto done;

    int64_t flow = overlap ? 0 : gf_dinic(N, src, snk, INT64_MAX);
    /* The flow is infinite iff it reaches the joint stand-in (see total). */
    GfCap vc = cap;
    if (split) vc.any_inf = cap.any_inf || vcap.any_inf;

    if (prop == GF_PROP_VALUE) {
        out = gf_value(&vc, flow);
    } else {
        /* Per-edge net flow, oriented. Edge k's forward arc sits at pos[k];
         * for an undirected edge without splitting, net flow a->b is
         * cap - res (negative: it runs b->a); with splitting the second
         * direction is the extra pair recorded in order. */
        int* ft = malloc((size_t)(m > 0 ? m : 1) * sizeof(int));
        int* fh = malloc((size_t)(m > 0 ? m : 1) * sizeof(int));
        int64_t* fv = malloc((size_t)(m > 0 ? m : 1) * sizeof(int64_t));
        long* fe = malloc((size_t)(m > 0 ? m : 1) * sizeof(long));
        if (!ft || !fh || !fv || !fe) { free(ft); free(fh); free(fv); free(fe); goto done; }
        long cnt = 0, p2 = m;
        for (long k = 0; k < m; k++) {
            int64_t f;
            if (edir[k] || !split) f = cap.c[k] - N->res[pos[k]];
            else { f = (cap.c[k] - N->res[pos[k]]) - (cap.c[k] - N->res[pos[p2]]); p2++; }
            if (f == 0 || overlap) continue;
            if (f > 0) { ft[cnt] = eu[k]; fh[cnt] = ev[k]; fv[cnt] = f; }
            else       { ft[cnt] = ev[k]; fh[cnt] = eu[k]; fv[cnt] = -f; }
            fe[cnt] = k; cnt++;
        }
        if (prop == GF_PROP_MATRIX) {
            out = gf_flow_matrix(n, ft, fh, fv, cnt, &vc);
        } else {
            /* Edges carrying flow, oriented along it, sorted by (tail, head)
             * -- the row-major order of the flow matrix, as Mathematica. */
            long* order = malloc((size_t)(cnt > 0 ? cnt : 1) * sizeof(long));
            long* bucket = calloc((size_t)n + 1, sizeof(long));
            long* tmp = malloc((size_t)(cnt > 0 ? cnt : 1) * sizeof(long));
            Expr** items = malloc((size_t)(cnt > 0 ? cnt : 1) * sizeof(Expr*));
            if (order && bucket && tmp && items) {
                /* LSD radix: by head, then stably by tail */
                for (long i = 0; i < cnt; i++) bucket[fh[i] + 1]++;
                for (int v = 0; v < n; v++) bucket[v + 1] += bucket[v];
                for (long i = 0; i < cnt; i++) tmp[bucket[fh[i]]++] = i;
                memset(bucket, 0, ((size_t)n + 1) * sizeof(long));
                for (long i = 0; i < cnt; i++) bucket[ft[i] + 1]++;
                for (int v = 0; v < n; v++) bucket[v + 1] += bucket[v];
                for (long i = 0; i < cnt; i++) order[bucket[ft[tmp[i]]]++] = tmp[i];
                Expr* hd = expr_new_symbol(SYM_DirectedEdge);
                Expr* hu = expr_new_symbol(SYM_UndirectedEdge);
                for (long i = 0; i < cnt; i++) {
                    long j = order[i];
                    items[i] = galg_make_edge(g, edir[fe[j]] ? hd : hu, ft[j], fh[j]);
                }
                expr_free(hd); expr_free(hu);
                out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)cnt);
            }
            free(order); free(bucket); free(tmp); free(items);
        }
        free(ft); free(fh); free(fv); free(fe);
    }

done:
    gf_net_free(N);
    free(pu); free(pv); free(cf); free(cb); free(pos);
    free(role); free(S); free(T);
    gf_cap_free(&cap); gf_cap_free(&vcap);
    return out;
}

/* ======================================================================== *
 * Global minimum cut
 * ======================================================================== */

/* Cut-family capacities: EdgeWeight when present, else 1. */
static int gf_weight_caps(const Expr* g, GfCap* cap) {
    if (g->data.function.arg_count == 3) {
        Expr* w = graph_resolve_edge_weights(g);
        if (!w) return 0;
        int ok = gf_caps_parse(w, galg_ne(g), cap);
        expr_free(w);
        return ok;
    }
    return gf_caps_parse(NULL, galg_ne(g), cap);
}

/* Undirected global min cut by Nagamochi-Ibaraki. Weighted edges
 * (eu[k], ev[k], w[k]) over n >= 2 vertices; parallel edges are summed. Writes
 * side[v] = 1 for one shore of a minimum cut and returns its value; -1 on
 * allocation failure. */
static int64_t gf_ni_mincut(int n, long m, const int* eu, const int* ev, const int64_t* w,
                            char* side) {
    int64_t result = -1;
    size_t nn = (size_t)n + 1, mm = (size_t)(2 * m) + 1;
    /* current contracted graph (CSR) */
    int* off = malloc(nn * sizeof(int));
    int* adj = malloc(mm * sizeof(int));
    int64_t* wt = malloc(mm * sizeof(int64_t));
    int* noff = malloc(nn * sizeof(int));
    int* nadj = malloc(mm * sizeof(int));
    int64_t* nwt = malloc(mm * sizeof(int64_t));
    int* owner = malloc(nn * sizeof(int));     /* original v -> current node */
    int* uf = malloc(nn * sizeof(int));
    int* newid = malloc(nn * sizeof(int));
    int64_t* key = malloc(nn * sizeof(int64_t));
    char* done = malloc(nn);
    int* order = malloc(nn * sizeof(int));
    int64_t* acc = malloc(nn * sizeof(int64_t));
    int* stamp = malloc(nn * sizeof(int));
    int* heapv = malloc((mm + nn) * sizeof(int));
    int64_t* heapk = malloc((mm + nn) * sizeof(int64_t));
    if (!off || !adj || !wt || !noff || !nadj || !nwt || !owner || !uf || !newid || !key
        || !done || !order || !acc || !stamp || !heapv || !heapk) goto out;

    /* initial graph: merge parallel edges by stamping */
    {
        int* cnt = calloc(nn, sizeof(int));
        int* tu = malloc(mm * sizeof(int));
        int64_t* tw = malloc(mm * sizeof(int64_t));
        if (!cnt || !tu || !tw) { free(cnt); free(tu); free(tw); goto out; }
        for (long k = 0; k < m; k++) if (eu[k] != ev[k]) { cnt[eu[k]]++; cnt[ev[k]]++; }
        int s = 0;
        for (int v = 0; v < n; v++) { off[v] = s; s += cnt[v]; }
        off[n] = s;
        for (int v = 0; v < n; v++) cnt[v] = off[v];
        for (long k = 0; k < m; k++) {
            int a = eu[k], b = ev[k];
            if (a == b) continue;
            tu[cnt[a]] = b; tw[cnt[a]++] = w[k];
            tu[cnt[b]] = a; tw[cnt[b]++] = w[k];
        }
        for (int v = 0; v < n; v++) stamp[v] = -1;
        int wpos = 0;
        for (int v = 0; v < n; v++) {
            int lo = off[v], hi = off[v + 1], start = wpos;
            off[v] = start;
            for (int j = lo; j < hi; j++) {
                int x = tu[j];
                if (stamp[x] >= start) wt[stamp[x]] += tw[j];
                else { stamp[x] = wpos; adj[wpos] = x; wt[wpos] = tw[j]; wpos++; }
            }
        }
        off[n] = wpos;
        free(cnt); free(tu); free(tw);
    }
    for (int v = 0; v < n; v++) owner[v] = v;

    int cn = n;                 /* current node count */
    int64_t lambda = INT64_MAX;
    long polls = 0;
    while (cn > 1) {
        if ((++polls & 7) == 0) galg_poll();
        /* bound by weighted degrees */
        int bestx = -1;
        for (int x = 0; x < cn; x++) {
            int64_t d = 0;
            for (int j = off[x]; j < off[x + 1]; j++) d += wt[j];
            if (d < lambda) { lambda = d; bestx = x; }
        }
        if (bestx >= 0) for (int v = 0; v < n; v++) side[v] = owner[v] == bestx;
        if (lambda == 0) break;

        /* maximum-adjacency order with a lazy max-heap; contract edges whose
         * attachment reaches lambda. */
        for (int x = 0; x < cn; x++) { key[x] = 0; done[x] = 0; uf[x] = x; }
        int hn = 0, cnt_order = 0;
        for (int start = 0; start < cn; start++) {
            if (done[start]) continue;
            /* push start (a new component of the current graph) */
            heapv[hn] = start; heapk[hn] = 0; hn++;
            while (hn > 0) {
                /* pop max */
                int x = heapv[0]; int64_t kx = heapk[0];
                hn--;
                if (hn > 0) {
                    int lv = heapv[hn]; int64_t lk = heapk[hn]; int i = 0;
                    for (;;) {
                        int c = 2 * i + 1;
                        if (c >= hn) break;
                        if (c + 1 < hn && heapk[c + 1] > heapk[c]) c++;
                        if (heapk[c] <= lk) break;
                        heapv[i] = heapv[c]; heapk[i] = heapk[c]; i = c;
                    }
                    heapv[i] = lv; heapk[i] = lk;
                }
                if (done[x] || kx != key[x]) continue;   /* stale */
                done[x] = 1;
                order[cnt_order++] = x;
                for (int j = off[x]; j < off[x + 1]; j++) {
                    int y = adj[j];
                    if (done[y]) continue;
                    key[y] += wt[j];
                    if (key[y] >= lambda) {
                        /* q(x,y) >= lambda: contract */
                        int a = x, b = y;
                        while (uf[a] != a) a = uf[a] = uf[uf[a]];
                        while (uf[b] != b) b = uf[b] = uf[uf[b]];
                        if (a != b) uf[a] = b;
                    }
                    /* push (y, key[y]) */
                    int i = hn++;
                    while (i > 0) {
                        int p = (i - 1) / 2;
                        if (heapk[p] >= key[y]) break;
                        heapv[i] = heapv[p]; heapk[i] = heapk[p]; i = p;
                    }
                    heapv[i] = y; heapk[i] = key[y];
                }
            }
        }
        /* the last two of the order: lambda(s, t) = its attachment = deg(t)
         * >= lambda, so they may always be merged (guarantees progress) */
        {
            int a = order[cn - 2], b = order[cn - 1];
            while (uf[a] != a) a = uf[a] = uf[uf[a]];
            while (uf[b] != b) b = uf[b] = uf[uf[b]];
            if (a != b) uf[a] = b;
        }
        /* relabel and rebuild */
        int nc = 0;
        for (int x = 0; x < cn; x++) {
            int r = x;
            while (uf[r] != r) r = uf[r];
            uf[x] = r;
        }
        for (int x = 0; x < cn; x++) newid[x] = -1;
        for (int x = 0; x < cn; x++) if (uf[x] == x) newid[x] = nc++;
        for (int x = 0; x < cn; x++) newid[x] = newid[uf[x]];
        /* group old nodes by new id */
        int* grp = order;           /* reuse: old nodes sorted by new id */
        {
            int* c2 = stamp;        /* reuse as counts */
            for (int y = 0; y <= nc; y++) c2[y] = 0;
            for (int x = 0; x < cn; x++) c2[newid[x] + 1]++;
            for (int y = 0; y < nc; y++) c2[y + 1] += c2[y];
            for (int x = 0; x < cn; x++) grp[c2[newid[x]]++] = x;
            /* c2[y] is now the end of group y; starts are c2[y-1] */
        }
        for (int y = 0; y < nc; y++) stamp[y] = -1;
        int wpos = 0, gi = 0;
        for (int y = 0; y < nc; y++) {
            int start = wpos;
            noff[y] = start;
            while (gi < cn && newid[grp[gi]] == y) {
                int x = grp[gi++];
                for (int j = off[x]; j < off[x + 1]; j++) {
                    int z = newid[adj[j]];
                    if (z == y) continue;
                    if (stamp[z] >= start) nwt[stamp[z]] += wt[j];
                    else { stamp[z] = wpos; nadj[wpos] = z; nwt[wpos] = wt[j]; wpos++; }
                }
            }
        }
        noff[nc] = wpos;
        { int* t = off; off = noff; noff = t; }
        { int* t = adj; adj = nadj; nadj = t; }
        { int64_t* t = wt; wt = nwt; nwt = t; }
        for (int v = 0; v < n; v++) owner[v] = newid[owner[v]];
        cn = nc;
    }
    result = lambda == INT64_MAX ? 0 : lambda;
out:
    free(off); free(adj); free(wt); free(noff); free(nadj); free(nwt); free(owner);
    free(uf); free(newid); free(key); free(done); free(order); free(acc); free(stamp);
    free(heapv); free(heapk);
    return result;
}

/* True iff every edge of g is undirected. */
static int gf_all_undirected(const Expr* g) { return graph_directed_edge_count(g) == 0; }

/* Directed global min cut: min over v != 0 of flow(0 -> v) and flow(v -> 0),
 * capped at the running bound. side[v] = 1 on the SOURCE side of the best
 * cut. Returns the value, -1 on failure. */
static int64_t gf_directed_mincut(const Expr* g, const GfCap* cap, char* side) {
    int n = galg_nv(g);
    GfNet* N = gf_edge_net(g, cap->c, NULL, 0, 0, NULL, NULL, NULL, 1);
    char* mark = malloc((size_t)n + 1);
    if (!N || !mark) { gf_net_free(N); free(mark); return -1; }
    int64_t best = INT64_MAX;
    for (int v = 1; v < n && best > 0; v++) {
        for (int dir = 0; dir < 2 && best > 0; dir++) {
            int s = dir ? v : 0, t = dir ? 0 : v;
            gf_net_reset(N);
            int64_t f = gf_dinic(N, s, t, best);
            if (f < best) {
                best = f;
                gf_reach_from(N, s, mark);
                memcpy(side, mark, (size_t)n);
            }
        }
    }
    gf_net_free(N); free(mark);
    return best == INT64_MAX ? 0 : best;
}

/* The minimum cut of g (n >= 2) with its capacities: value in *val, side[]
 * one shore (the source side when directed). Returns 1, 0 on failure. */
static int gf_global_mincut(const Expr* g, const GfCap* cap, int64_t* val, char* side) {
    int n = galg_nv(g);
    long m = galg_ne(g);
    if (gf_all_undirected(g)) {
        const int *eu, *ev;
        const unsigned char* edir;
        graph_edge_indices(g, &eu, &ev, &edir);
        *val = gf_ni_mincut(n, m, eu, ev, cap->c, side);
    } else {
        *val = gf_directed_mincut(g, cap, side);
    }
    return *val >= 0;
}

Expr* builtin_find_minimum_cut(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    int n = galg_nv(g);
    if (n < 2) return NULL;                         /* Mathematica: unevaluated */
    GfCap cap;
    if (!gf_weight_caps(g, &cap)) return NULL;
    char* side = calloc((size_t)n, 1);
    int* a = malloc((size_t)n * sizeof(int));
    int* b = malloc((size_t)n * sizeof(int));
    Expr* out = NULL;
    int64_t val;
    if (side && a && b && gf_global_mincut(g, &cap, &val, side)) {
        int directed = !gf_all_undirected(g);
        /* first listed part: the source side (directed) or the side without
         * vertex 0 (undirected) */
        char want = directed ? 1 : (char)!side[0];
        int na = 0, nb = 0;
        for (int v = 0; v < n; v++) { if (side[v] == want) a[na++] = v; else b[nb++] = v; }
        Expr* parts[2] = { galg_vertex_list(g, a, na, 0), galg_vertex_list(g, b, nb, 0) };
        Expr* items[2] = { gf_value(&cap, val),
                           expr_new_function(expr_new_symbol(SYM_List), parts, 2) };
        out = expr_new_function(expr_new_symbol(SYM_List), items, 2);
    }
    free(side); free(a); free(b); gf_cap_free(&cap);
    return out;
}

/* Edges of g crossing from side 1 to side 0 (undirected edges crossing either
 * way), in EdgeList order. */
static Expr* gf_cut_edges(const Expr* g, const char* side) {
    const int *eu, *ev;
    const unsigned char* edir;
    graph_edge_indices(g, &eu, &ev, &edir);
    long m = galg_ne(g), cnt = 0;
    long* idx = malloc((size_t)(m > 0 ? m : 1) * sizeof(long));
    if (!idx) return NULL;
    for (long k = 0; k < m; k++) {
        int a = side[eu[k]], b = side[ev[k]];
        if (edir[k] ? (a && !b) : (a != b)) idx[cnt++] = k;
    }
    Expr* out = galg_edge_list(g, idx, cnt);
    free(idx);
    return out;
}

/* s-t edge cut with g's cut capacities: side[] = residual reach from s
 * (the cut closest to s). Returns the flow value, -1 on failure. */
static int64_t gf_st_cut(const Expr* g, const GfCap* cap, int s, int t, char* side) {
    GfNet* N = gf_edge_net(g, cap->c, NULL, 0, 0, NULL, NULL, NULL, 0);
    if (!N) return -1;
    int64_t f = gf_dinic(N, s, t, INT64_MAX);
    if (side) gf_reach_from(N, s, side);
    gf_net_free(N);
    return f;
}

Expr* builtin_find_edge_cut(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 3) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    int n = galg_nv(g);
    GfCap cap;
    if (!gf_weight_caps(g, &cap)) return NULL;
    char* side = calloc((size_t)(n > 0 ? n : 1), 1);
    Expr* out = NULL;
    if (!side) goto done;
    if (argc == 1) {
        if (n < 2) { out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0); goto done; }
        int64_t val;
        if (!gf_global_mincut(g, &cap, &val, side)) goto done;
        out = gf_cut_edges(g, side);
    } else {
        int s = galg_vertex_arg(g, res->data.function.args[1]);
        int t = galg_vertex_arg(g, res->data.function.args[2]);
        if (s < 0 || t < 0 || s == t) goto done;
        if (gf_st_cut(g, &cap, s, t, side) < 0) goto done;
        out = gf_cut_edges(g, side);
    }
done:
    free(side); gf_cap_free(&cap);
    return out;
}

Expr* builtin_edge_connectivity(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 3) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    int n = galg_nv(g);
    GfCap cap;
    if (!gf_weight_caps(g, &cap)) return NULL;
    Expr* out = NULL;
    char* side = calloc((size_t)(n > 0 ? n : 1), 1);
    if (!side) goto done;
    if (argc == 1) {
        if (n < 2) goto done;                       /* Mathematica: unevaluated */
        int64_t val;
        if (gf_global_mincut(g, &cap, &val, side)) out = gf_value(&cap, val);
    } else {
        int s = galg_vertex_arg(g, res->data.function.args[1]);
        int t = galg_vertex_arg(g, res->data.function.args[2]);
        if (s < 0 || t < 0 || s == t) goto done;
        int64_t f = gf_st_cut(g, &cap, s, t, NULL);
        if (f >= 0) out = gf_value(&cap, f);
    }
done:
    free(side); gf_cap_free(&cap);
    return out;
}

/* ======================================================================== *
 * Vertex connectivity and separators (underlying undirected graph)
 * ======================================================================== */

/* Split network of the simple undirected graph u: v_in = 2v, v_out = 2v+1,
 * v_in -> v_out capacity 1, each edge {a,b} as a_out -> b_in and
 * b_out -> a_in with capacity n (effectively infinite). */
static GfNet* gf_split_net(const GalgUG* u) {
    int n = u->n;
    long np = n + 2 * u->m;
    int* pu = malloc((size_t)(np > 0 ? np : 1) * sizeof(int));
    int* pv = malloc((size_t)(np > 0 ? np : 1) * sizeof(int));
    int64_t* cf = malloc((size_t)(np > 0 ? np : 1) * sizeof(int64_t));
    int64_t* cb = calloc((size_t)(np > 0 ? np : 1), sizeof(int64_t));
    GfNet* N = NULL;
    if (pu && pv && cf && cb) {
        long p = 0;
        for (int v = 0; v < n; v++) { pu[p] = 2 * v; pv[p] = 2 * v + 1; cf[p] = 1; p++; }
        for (int a = 0; a < n; a++)
            for (int j = u->off[a]; j < u->off[a + 1]; j++) {
                int b = u->adj[j];
                pu[p] = 2 * a + 1; pv[p] = 2 * b; cf[p] = n; p++;
            }
        N = gf_net_build(2 * n, p, pu, pv, cf, cb, NULL, 1);
    }
    free(pu); free(pv); free(cf); free(cb);
    return N;
}

static int gf_ug_adjacent(const GalgUG* u, int a, int b) {
    int lo = u->off[a], hi = u->off[a + 1] - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2, x = u->adj[mid];
        if (x == b) return 1;
        if (x < b) lo = mid + 1; else hi = mid - 1;
    }
    return 0;
}

/* Local vertex connectivity kappa(s, t) (s, t non-adjacent), capped at limit.
 * Leaves the flow in N. */
static int64_t gf_local_kappa(GfNet* N, int s, int t, int64_t limit) {
    gf_net_reset(N);
    return gf_dinic(N, 2 * s + 1, 2 * t, limit);
}

/* Separator from the flow in N: vertices whose in-copy cannot reach t_in but
 * whose out-copy can (the cut closest to t). Writes into idx, returns count. */
static int gf_separator(GfNet* N, int n, int t, int* idx, char* mark) {
    gf_reach_to(N, 2 * t, mark);
    int k = 0;
    for (int v = 0; v < n; v++) if (!mark[2 * v] && mark[2 * v + 1]) idx[k++] = v;
    return k;
}

/* Global vertex connectivity (Esfahanian-Hakimi). If sep != NULL, writes a
 * minimum separator there (count in *nsep). u must have n >= 1. */
static long gf_kappa_global(const GalgUG* u, int* sep, int* nsep) {
    int n = u->n;
    if (nsep) *nsep = 0;
    if (n <= 1) return 0;
    /* complete graph: n - 1 (separator: the first n-1 vertices) */
    if (u->m == (long)n * (n - 1) / 2) {
        if (sep) { for (int v = 0; v < n - 1; v++) sep[v] = v; *nsep = n - 1; }
        return n - 1;
    }
    GfNet* N = gf_split_net(u);
    char* mark = malloc((size_t)(2 * n) + 1);
    if (!N || !mark) { gf_net_free(N); free(mark); return -1; }
    int v0 = 0;
    for (int v = 1; v < n; v++) if (u->off[v + 1] - u->off[v] < u->off[v0 + 1] - u->off[v0]) v0 = v;
    int64_t best = n - 1;
    int bs = -1, bt = -1;
    char* nb = calloc((size_t)n, 1);
    if (!nb) { gf_net_free(N); free(mark); return -1; }
    for (int j = u->off[v0]; j < u->off[v0 + 1]; j++) nb[u->adj[j]] = 1;
    for (int w = 0; w < n && best > 0; w++) {
        if (w == v0 || nb[w]) continue;
        int64_t k = gf_local_kappa(N, v0, w, best);
        if (k < best) { best = k; bs = v0; bt = w; }
    }
    for (int i = u->off[v0]; i < u->off[v0 + 1] && best > 0; i++)
        for (int j = i + 1; j < u->off[v0 + 1] && best > 0; j++) {
            int x = u->adj[i], y = u->adj[j];
            if (gf_ug_adjacent(u, x, y)) continue;
            int64_t k = gf_local_kappa(N, x, y, best);
            if (k < best) { best = k; bs = x; bt = y; }
        }
    if (sep) {
        if (bs < 0) {
            /* no pair beat n-1: only possible for complete graphs (handled) */
            for (int v = 0; v < n - 1; v++) sep[v] = v;
            *nsep = n - 1;
        } else if (best == 0) {
            *nsep = 0;
        } else {
            gf_local_kappa(N, bs, bt, INT64_MAX);
            *nsep = gf_separator(N, n, bt, sep, mark);
        }
    }
    free(nb); free(mark); gf_net_free(N);
    return best;
}

long galg_vertex_connectivity(const Expr* g, int s, int t) {
    GalgUG* u = galg_ug_from_graph(g);
    if (!u) return -1;
    long r;
    if (s < 0 || t < 0) {
        r = gf_kappa_global(u, NULL, NULL);
    } else if (s == t || gf_ug_adjacent(u, s, t)) {
        r = 0;
    } else {
        GfNet* N = gf_split_net(u);
        r = N ? gf_local_kappa(N, s, t, INT64_MAX) : -1;
        gf_net_free(N);
    }
    galg_ug_free(u);
    return r;
}

Expr* builtin_find_vertex_cut(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 3) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    GalgUG* u = galg_ug_from_graph(g);
    if (!u) return NULL;
    int n = u->n;
    int* sep = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    Expr* out = NULL;
    if (!sep) goto done;
    if (argc == 1) {
        int ns = 0;
        /* Mathematica: a graph with a directed edge whose underlying graph is
         * complete (n >= 3) has no vertex cut at all -> {} (an undirected
         * complete graph gets its first n-1 vertices). */
        if (n >= 3 && u->m == (long)n * (n - 1) / 2 && graph_directed_edge_count(g) > 0) {
            out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
            goto done;
        }
        if (gf_kappa_global(u, sep, &ns) < 0) goto done;
        out = galg_vertex_list(g, sep, ns, 1);
    } else {
        int s = galg_vertex_arg(g, res->data.function.args[1]);
        int t = galg_vertex_arg(g, res->data.function.args[2]);
        if (s < 0 || t < 0 || s == t) goto done;
        if (gf_ug_adjacent(u, s, t)) { out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0); goto done; }
        GfNet* N = gf_split_net(u);
        char* mark = malloc((size_t)(2 * n) + 1);
        if (N && mark) {
            gf_local_kappa(N, s, t, INT64_MAX);
            int ns = gf_separator(N, n, t, sep, mark);
            out = galg_vertex_list(g, sep, ns, 1);
        }
        free(mark); gf_net_free(N);
    }
done:
    free(sep); galg_ug_free(u);
    return out;
}

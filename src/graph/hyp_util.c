/* hyp_util.c - Hypergraph construction, validation, the validated-hypergraph
 * memo, HypergraphQ, and the accessor dispatch shared with the Graph heads
 * (VertexList, EdgeList, VertexCount, EdgeCount, VertexDegree,
 * IncidenceMatrix).
 *
 * See graph_hyper.h for the representation and the ordered/unordered choice.
 *
 * CONSTRUCTION. Accepted forms:
 *   Hypergraph[{e1, e2, ...}]          vertices derived in first-appearance order
 *   Hypergraph[{v1, ...}, {e1, ...}]   explicit vertices (duplicates dropped,
 *                                      first occurrence kept); every hyperedge
 *                                      element must be one of them
 *   Hypergraph[g]                      a Graph: each edge u<->v / u->v becomes
 *                                      the hyperedge {u, v}
 * Every hyperedge must be a List. Anything else is left unevaluated. A valid
 * hypergraph is its own fixed point (the constructor returns NULL for it).
 *
 * THE MEMO. graph_util.c's design, re-implemented locally (graph.h's memo is
 * Graph-specific): the last HYP_MEMO_SLOTS
 * valid hypergraphs, keyed by node pointer, each slot holding a reference so
 * the node stays alive -- hence unrecyclable and immutable (mutators unshare a
 * node with refcount > 1). A slot stores the vertex index, hyperedges as a
 * vertex-index CSR (raw and de-duplicated), and, lazily, the
 * vertex->hyperedge incidence CSR. Invalid objects are never memoized.
 *
 * Memory (SPEC section 4): builtins take `res`, return a fresh result or NULL.
 */

#include "graph_hyper.h"
#include "graph.h"
#include "expr.h"
#include "pack.h"
#include "sym_names.h"
#include "sym_intern.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ---- small helpers -------------------------------------------------------- */

static const char* s_hypergraph = NULL;
const char* hyp_sym_hypergraph(void) {
    if (!s_hypergraph) s_hypergraph = intern_symbol("Hypergraph");
    return s_hypergraph;
}
#define sym_hypergraph hyp_sym_hypergraph

static int head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

/* Hypergraph[List, List] shape (not yet validated). */
static int hyp_shape_ok(const Expr* h) {
    if (!head_is(h, sym_hypergraph()) || h->data.function.arg_count != 2) return 0;
    return graph_is_list(h->data.function.args[0]) && graph_is_list(h->data.function.args[1]);
}

/* ---- the memo ------------------------------------------------------------- */
#define HYP_MEMO_SLOTS 4

typedef struct {
    Expr*      h;        /* owned reference; NULL = empty slot              */
    GraphVIdx* ix;       /* vertex -> index; keys borrowed from h's verts    */
    int n, m;
    int *eoff, *ev;      /* raw hyperedges as vertex indices                 */
    int *soff, *sv;      /* distinct-vertex hyperedges (alias raw if no rep) */
    int has_rep;
    int *voff, *ve;      /* incidence, lazily built; NULL until then         */
    uint64_t last_used;  /* LRU clock of the latest hit/insert               */
} HypMemo;

/* Least-recently-used eviction, like the graph memo (graph_util.c), so a
 * cyclic pattern over a few live hypergraphs does not thrash. Module-static
 * with no locking, following the same g_qqbar_cache precedent: hypergraph
 * builtins run only on the evaluator thread (gmet worker threads never call
 * into Expr code). */
static HypMemo g_hyp_memo[HYP_MEMO_SLOTS];
static uint64_t g_hyp_memo_clock = 0;

static void hyp_memo_clear(HypMemo* s) {
    graph_vidx_free(s->ix);
    if (s->has_rep) { free(s->soff); free(s->sv); }
    free(s->eoff); free(s->ev);
    free(s->voff); free(s->ve);
    if (s->h) expr_free(s->h);
    memset(s, 0, sizeof(*s));
}

/* ---- integer fast path ----------------------------------------------------
 * Most hypergraphs have machine-integer vertices (Range[n], Wolfram-model
 * states). Resolving each of the sum|e| hyperedge elements through the
 * expr_hash index was the dominant cost of construction, so when every vertex
 * is an EXPR_INTEGER in a range at most 4n + 1024 wide, elements resolve by
 * direct addressing instead. Exact: a canonical expression SameQ to a machine
 * integer is that EXPR_INTEGER (small BigInts are normalised), so an element
 * of any other type is simply not a vertex. The hash index is still built over
 * the n vertices for single lookups (VertexDegree[h, v], ...). */
typedef struct { int64_t lo; size_t span; int* at; } IntMap;

/* Map from the n vertex nodes vs[0..n-1] (value -> position, first wins), or
 * at == NULL when the fast path does not apply. */
static IntMap intmap_build(Expr* const* vs, size_t n) {
    IntMap im = { 0, 0, NULL };
    if (n == 0) return im;
    int64_t lo = INT64_MAX, hi = INT64_MIN;
    for (size_t i = 0; i < n; i++) {
        if (vs[i]->type != EXPR_INTEGER) return im;
        int64_t x = vs[i]->data.integer;
        if (x < lo) lo = x;
        if (x > hi) hi = x;
    }
    /* Span in unsigned arithmetic: hi - lo as int64_t overflows (UB) when the
     * vertices straddle more than INT64_MAX. Once it is bounded below, every
     * x - lo for x in [lo, hi] is small, so the signed subtractions are safe. */
    uint64_t width = (uint64_t)hi - (uint64_t)lo;
    if (width > (uint64_t)n * 4 + 1024) return im;
    size_t span = (size_t)width + 1;
    int* at = malloc(span * sizeof(int));
    if (!at) return im;
    for (size_t k = 0; k < span; k++) at[k] = -1;
    for (size_t i = 0; i < n; i++) {
        size_t k = (size_t)(vs[i]->data.integer - lo);
        if (at[k] < 0) at[k] = (int)i;
    }
    im.lo = lo; im.span = span; im.at = at;
    return im;
}

static inline int intmap_get(const IntMap* im, const Expr* x) {
    if (x->type != EXPR_INTEGER) return -1;
    int64_t d = x->data.integer - im->lo;
    if (d < 0 || (uint64_t)d >= im->span) return -1;
    return im->at[d];
}

/* Given an index over h's (distinct) vertex list, resolve every hyperedge and
 * store the memo entry. Takes ownership of ix. NULL if a hyperedge is not a
 * List or names a non-vertex. */
static HypMemo* hyp_memo_from(const Expr* h, GraphVIdx* ix) {
    const Expr* verts = h->data.function.args[0];
    const Expr* edges = h->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t m = edges->data.function.arg_count;
    size_t tot = 0;
    for (size_t j = 0; j < m; j++) {
        const Expr* e = edges->data.function.args[j];
        if (!graph_is_list(e)) { graph_vidx_free(ix); return NULL; }
        tot += e->data.function.arg_count;
    }
    if (tot > (size_t)INT32_MAX || m > (size_t)INT32_MAX - 1) { graph_vidx_free(ix); return NULL; }
    int* eoff = malloc((m + 1) * sizeof(int));
    int* ev = malloc((tot ? tot : 1) * sizeof(int));
    int* stamp = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!eoff || !ev || !stamp) goto fail;
    for (int i = 0; i < n; i++) stamp[i] = -1;
    IntMap im = intmap_build(verts->data.function.args, (size_t)n);
    int has_rep = 0;
    size_t p = 0;
    for (size_t j = 0; j < m; j++) {
        const Expr* e = edges->data.function.args[j];
        eoff[j] = (int)p;
        for (size_t k = 0; k < e->data.function.arg_count; k++) {
            int vi = im.at ? intmap_get(&im, e->data.function.args[k])
                           : graph_vidx_get(ix, e->data.function.args[k]);
            if (vi < 0) { free(im.at); goto fail; }
            if (stamp[vi] == (int)j) has_rep = 1;
            stamp[vi] = (int)j;
            ev[p++] = vi;
        }
    }
    eoff[m] = (int)p;
    free(im.at);

    int *soff = eoff, *sv = ev;
    if (has_rep) {
        soff = malloc((m + 1) * sizeof(int));
        sv = malloc((tot ? tot : 1) * sizeof(int));
        if (!soff || !sv) { free(soff); free(sv); goto fail; }
        for (int i = 0; i < n; i++) stamp[i] = -1;
        size_t q = 0;
        for (size_t j = 0; j < m; j++) {
            soff[j] = (int)q;
            for (int k = eoff[j]; k < eoff[j + 1]; k++) {
                int vi = ev[k];
                if (stamp[vi] == (int)j) continue;
                stamp[vi] = (int)j;
                sv[q++] = vi;
            }
        }
        soff[m] = (int)q;
    }
    free(stamp);

    HypMemo* s = &g_hyp_memo[0];
    for (int i = 0; i < HYP_MEMO_SLOTS; i++) {
        if (!g_hyp_memo[i].h) { s = &g_hyp_memo[i]; break; }
        if (g_hyp_memo[i].last_used < s->last_used) s = &g_hyp_memo[i];
    }
    if (s->h) hyp_memo_clear(s);
    s->h = expr_copy((Expr*)h);        /* refcount bump only; never mutates */
    s->ix = ix;
    s->n = n; s->m = (int)m;
    s->eoff = eoff; s->ev = ev;
    s->soff = soff; s->sv = sv;
    s->has_rep = has_rep;
    s->voff = NULL; s->ve = NULL;
    s->last_used = ++g_hyp_memo_clock;
    return s;

fail:
    free(eoff); free(ev); free(stamp);
    graph_vidx_free(ix);
    return NULL;
}

/* The memo entry for h, validating (and memoizing) on a miss; NULL if invalid. */
static HypMemo* hyp_memo(const Expr* h) {
    for (int i = 0; i < HYP_MEMO_SLOTS; i++)
        if (g_hyp_memo[i].h == h) {
            g_hyp_memo[i].last_used = ++g_hyp_memo_clock;
            return &g_hyp_memo[i];
        }
    if (!hyp_shape_ok(h)) return NULL;
    const Expr* verts = h->data.function.args[0];
    size_t n = verts->data.function.arg_count;
    if (n > (size_t)INT32_MAX) return NULL;
    GraphVIdx* ix = graph_vidx_new(n);
    if (!ix) return NULL;
    for (size_t i = 0; i < n; i++)
        if (!graph_vidx_put(ix, verts->data.function.args[i], (int)i)) {
            graph_vidx_free(ix);            /* repeated vertex: not canonical */
            return NULL;
        }
    return hyp_memo_from(h, ix);
}

/* Build the vertex -> hyperedge incidence CSR (from the distinct sets, so each
 * hyperedge is listed once per vertex, in ascending hyperedge order). */
static int hyp_memo_incidence(HypMemo* s) {
    if (s->voff) return 1;
    int n = s->n, m = s->m;
    int tot = s->soff[m];
    int* voff = calloc((size_t)n + 1, sizeof(int));
    int* ve = malloc((size_t)(tot ? tot : 1) * sizeof(int));
    if (!voff || !ve) { free(voff); free(ve); return 0; }
    for (int k = 0; k < tot; k++) voff[s->sv[k] + 1]++;
    for (int i = 0; i < n; i++) voff[i + 1] += voff[i];
    int* fill = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!fill) { free(voff); free(ve); return 0; }
    memcpy(fill, voff, (size_t)n * sizeof(int));
    for (int j = 0; j < m; j++)
        for (int k = s->soff[j]; k < s->soff[j + 1]; k++)
            ve[fill[s->sv[k]]++] = j;
    free(fill);
    s->voff = voff; s->ve = ve;
    return 1;
}

int hypergraph_is_valid(const Expr* h) {
    if (!head_is(h, sym_hypergraph())) return 0;
    return hyp_memo(h) != NULL;
}

int hyp_view(const Expr* h, HypView* v, int want_incidence) {
    if (!head_is(h, sym_hypergraph())) return 0;
    HypMemo* s = hyp_memo(h);
    if (!s) return 0;
    if (want_incidence && !hyp_memo_incidence(s)) return 0;
    v->n = s->n; v->m = s->m;
    v->verts = h->data.function.args[0];
    v->edges = h->data.function.args[1];
    v->eoff = s->eoff; v->ev = s->ev;
    v->soff = s->soff; v->sv = s->sv;
    v->voff = s->voff; v->ve = s->ve;
    v->ix = s->ix;
    return 1;
}

/* ---- construction --------------------------------------------------------- */

static Expr* list_of(Expr** items, size_t n) {
    return expr_new_function(expr_new_symbol(SYM_List), items, n);
}

/* Hypergraph[g] for a valid Graph g. */
static Expr* hyp_from_graph(const Expr* g) {
    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    size_t m = edges->data.function.arg_count;
    Expr** es = calloc(m ? m : 1, sizeof(Expr*));
    if (!es) return NULL;
    Expr* lh = expr_new_symbol(SYM_List);
    for (size_t j = 0; j < m; j++) {
        const Expr* e = edges->data.function.args[j];
        Expr* a[2] = { expr_copy(e->data.function.args[0]), expr_copy(e->data.function.args[1]) };
        es[j] = expr_new_function(expr_copy(lh), a, 2);
    }
    expr_free(lh);
    Expr* args[2] = { expr_copy((Expr*)verts), list_of(es, m) };
    free(es);
    return expr_new_function(expr_new_symbol(sym_hypergraph()), args, 2);
}

Expr* builtin_hypergraph(Expr* res) {
    if (hypergraph_is_valid(res)) return NULL;      /* canonical: fixed point */
    size_t argc = res->data.function.arg_count;
    if (argc == 1) {
        const Expr* a = res->data.function.args[0];
        if (graph_is_valid(a)) return hyp_from_graph(a);
        if (!graph_is_list(a)) return NULL;
        /* Derive vertices in first-appearance order. */
        size_t m = a->data.function.arg_count, tot = 0;
        for (size_t j = 0; j < m; j++) {
            const Expr* e = a->data.function.args[j];
            if (!graph_is_list(e)) return NULL;
            tot += e->data.function.arg_count;
        }
        if (tot > (size_t)INT32_MAX) return NULL;
        Expr** vs = malloc((tot ? tot : 1) * sizeof(Expr*));
        if (!vs) return NULL;
        size_t nv = 0;
        GraphVIdx* ix = NULL;
        /* Integer fast path (see IntMap): first appearance by direct
         * addressing, then the hash index over the distinct vertices only. */
        int64_t lo = INT64_MAX, hi = INT64_MIN;
        int allint = tot > 0;
        for (size_t j = 0; j < m && allint; j++) {
            const Expr* e = a->data.function.args[j];
            for (size_t k = 0; k < e->data.function.arg_count; k++) {
                const Expr* x = e->data.function.args[k];
                if (x->type != EXPR_INTEGER) { allint = 0; break; }
                if (x->data.integer < lo) lo = x->data.integer;
                if (x->data.integer > hi) hi = x->data.integer;
            }
        }
        unsigned char* seen = NULL;
        if (allint && m > 0 && hi >= lo                      /* unsigned span: no UB */
            && (uint64_t)hi - (uint64_t)lo <= (uint64_t)tot * 4 + 1024)
            seen = calloc((size_t)((uint64_t)hi - (uint64_t)lo) + 1, 1);
        if (seen) {
            for (size_t j = 0; j < m; j++) {
                const Expr* e = a->data.function.args[j];
                for (size_t k = 0; k < e->data.function.arg_count; k++) {
                    Expr* x = e->data.function.args[k];
                    size_t d = (size_t)(x->data.integer - lo);
                    if (!seen[d]) { seen[d] = 1; vs[nv++] = expr_copy(x); }
                }
            }
            free(seen);
            ix = graph_vidx_new(nv);
            if (ix) for (size_t i = 0; i < nv; i++) graph_vidx_put(ix, vs[i], (int)i);
        } else {
            ix = graph_vidx_new(tot);
            if (ix)
                for (size_t j = 0; j < m; j++) {
                    const Expr* e = a->data.function.args[j];
                    for (size_t k = 0; k < e->data.function.arg_count; k++) {
                        Expr* x = e->data.function.args[k];
                        if (graph_vidx_put(ix, x, (int)nv)) vs[nv++] = expr_copy(x);
                    }
                }
        }
        if (!ix) {
            for (size_t i = 0; i < nv; i++) expr_free(vs[i]);
            free(vs);
            return NULL;
        }
        Expr* args[2] = { list_of(vs, nv), expr_copy((Expr*)a) };
        free(vs);
        Expr* h = expr_new_function(expr_new_symbol(sym_hypergraph()), args, 2);
        /* ix's keys are the edge-element nodes, which h keeps alive (the same
         * pointers sit in both its vertex List and its edge List). */
        if (!hyp_memo_from(h, ix)) { expr_free(h); return NULL; }
        return h;
    }
    if (argc != 2) return NULL;
    const Expr* vin = res->data.function.args[0];
    const Expr* ein = res->data.function.args[1];
    if (!graph_is_list(vin) || !graph_is_list(ein)) return NULL;
    size_t n = vin->data.function.arg_count;
    if (n > (size_t)INT32_MAX) return NULL;
    GraphVIdx* ix = graph_vidx_new(n);
    Expr** vs = malloc((n ? n : 1) * sizeof(Expr*));
    if (!ix || !vs) { graph_vidx_free(ix); free(vs); return NULL; }
    size_t nv = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* x = vin->data.function.args[i];
        if (graph_vidx_put(ix, x, (int)nv)) vs[nv++] = expr_copy(x);
    }
    if (nv == n) {
        /* No duplicate vertices, so res itself failed on its hyperedges:
         * nothing canonical to build. */
        for (size_t i = 0; i < nv; i++) expr_free(vs[i]);
        free(vs); graph_vidx_free(ix);
        return NULL;
    }
    Expr* args[2] = { list_of(vs, nv), expr_copy((Expr*)ein) };
    free(vs);
    Expr* h = expr_new_function(expr_new_symbol(sym_hypergraph()), args, 2);
    if (!hyp_memo_from(h, ix)) { expr_free(h); return NULL; }
    return h;
}

Expr* builtin_hypergraph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    return expr_new_symbol(hypergraph_is_valid(res->data.function.args[0]) ? SYM_True : SYM_False);
}

/* ---- accessor dispatch ---------------------------------------------------- */

Expr* hyp_vertex_list(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* h = res->data.function.args[0];
    if (!hypergraph_is_valid(h)) return NULL;
    return expr_copy(h->data.function.args[0]);
}

Expr* hyp_edge_list(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* h = res->data.function.args[0];
    if (!hypergraph_is_valid(h)) return NULL;
    return expr_copy(h->data.function.args[1]);
}

Expr* hyp_vertex_count(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* h = res->data.function.args[0];
    if (!hypergraph_is_valid(h)) return NULL;
    return expr_new_integer((int64_t)h->data.function.args[0]->data.function.arg_count);
}

Expr* hyp_edge_count(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* h = res->data.function.args[0];
    if (!hypergraph_is_valid(h)) return NULL;
    return expr_new_integer((int64_t)h->data.function.args[1]->data.function.arg_count);
}

/* A List of machine integers: packed (int64 buffer) above the packing
 * threshold, else a plain List. */
Expr* hyp_int_list(const int64_t* xs, size_t n) {
    int64_t* buf = NULL;
    Expr* packed = ndbuild_open_i64((int64_t)n, &buf);
    if (packed) {
        if (n) memcpy(buf, xs, n * sizeof(int64_t));
        return packed;
    }
    Expr** items = malloc((n ? n : 1) * sizeof(Expr*));
    if (!items) return NULL;
    for (size_t i = 0; i < n; i++) items[i] = expr_new_integer(xs[i]);
    Expr* out = list_of(items, n);
    free(items);
    return out;
}

/* VertexDegree[h]: for each vertex, the number of hyperedges containing it
 * (repeated hyperedges count separately; a vertex repeated inside one
 * hyperedge counts once). VertexDegree[h, v]: that number for v. */
Expr* hyp_vertex_degree(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 1)) return NULL;
    if (argc == 2) {
        int i = graph_vidx_get(v.ix, res->data.function.args[1]);
        if (i < 0) return NULL;
        return expr_new_integer((int64_t)(v.voff[i + 1] - v.voff[i]));
    }
    int64_t* d = malloc((size_t)(v.n > 0 ? v.n : 1) * sizeof(int64_t));
    if (!d) return NULL;
    for (int i = 0; i < v.n; i++) d[i] = v.voff[i + 1] - v.voff[i];
    Expr* out = hyp_int_list(d, (size_t)v.n);
    free(d);
    return out;
}

/* IncidenceMatrix[h]: the n x m matrix whose (i, j) entry is the number of
 * times vertex i occurs in hyperedge j -- 0/1 for set-like hyperedges, and the
 * multiplicity otherwise, as IncidenceMatrix gives 2 for a Graph self-loop.
 * Packed int64 above the packing threshold. */
Expr* hyp_incidence_matrix(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 0)) return NULL;
    size_t n = (size_t)v.n, m = (size_t)v.m;
    if (m && n > SIZE_MAX / 8 / m) return NULL;
    int64_t dims[2] = { (int64_t)n, (int64_t)m };
    void* raw = NULL;
    Expr* packed = (n && m) ? ndbuild_open(2, dims, NDT_INT64, &raw) : NULL;
    if (packed) {
        int64_t* buf = raw;
        memset(buf, 0, n * m * sizeof(int64_t));
        for (size_t j = 0; j < m; j++)
            for (int k = v.eoff[j]; k < v.eoff[j + 1]; k++)
                buf[(size_t)v.ev[k] * m + j]++;
        return packed;
    }
    int* grid = calloc((n * m) > 0 ? n * m : 1, sizeof(int));
    Expr** rows = malloc((n ? n : 1) * sizeof(Expr*));
    Expr** row = malloc((m ? m : 1) * sizeof(Expr*));
    if (!grid || !rows || !row) { free(grid); free(rows); free(row); return NULL; }
    for (size_t j = 0; j < m; j++)
        for (int k = v.eoff[j]; k < v.eoff[j + 1]; k++)
            grid[(size_t)v.ev[k] * m + j]++;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) row[j] = expr_new_integer(grid[i * m + j]);
        rows[i] = list_of(row, m);
    }
    Expr* out = list_of(rows, n);
    free(grid); free(rows); free(row);
    return out;
}

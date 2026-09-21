/* graph_util.c - shared, read-only helpers for the graph subsystem.
 *
 * Graphs are ordinary Expr trees; these helpers inspect the canonical form
 *
 *     Graph[ List[v1, ...], List[edge1, ...] ]
 *
 * where each edge is a 2-argument DirectedEdge[u, v] or UndirectedEdge[u, v].
 * Nothing here allocates or mutates; ownership contracts live in the callers.
 *
 * Vertex membership resolves through an expr_hash index (GraphVIdx below), built once
 * per validation or adjacency pass. The MVP did a linear expr_eq scan per
 * lookup, which made validation O(E*V) and parallel-edge detection O(E^2): a
 * 20000-vertex, 40000-edge graph spent ~14 s inside Graph[], and paid it again
 * in every accessor, since they all begin with graph_is_valid. Both passes are
 * now O(V + E) expected.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

/* True iff e is a function node whose head is the interned symbol `sym`. */
static int head_is_sym(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

int graph_is_list(const Expr* e) {
    return head_is_sym(e, SYM_List);
}

const char* graph_edge_kind(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2)
        return NULL;
    if (head_is_sym(e, SYM_DirectedEdge))   return SYM_DirectedEdge;
    if (head_is_sym(e, SYM_UndirectedEdge)) return SYM_UndirectedEdge;
    return NULL;
}

/* ---- Vertex index --------------------------------------------------------- *
 * Maps a vertex Expr to its position in the canonical vertex List. Vertices are
 * arbitrary expressions, so equality is expr_eq and the hash must therefore be
 * expr_hash -- the two agree by contract (see the "identity trio" note in
 * src/expr.c). Open addressing with linear probing, kept at load factor <= 0.5.
 *
 * Keys are borrowed pointers into the graph's vertex List; the index is only
 * ever used within a single call, while that graph is alive.
 *
 * A repeated vertex keeps its FIRST index, matching the linear scan in
 * graph_vertex_index that this replaces. */
struct GraphVIdx {
    const Expr** key;   /* NULL slot = empty; keys are borrowed               */
    int*         idx;
    size_t       mask;  /* capacity - 1; capacity is a power of two           */
    size_t       count;
};

void graph_vidx_free(GraphVIdx* ix) {
    if (!ix) return;
    free(ix->key);
    free(ix->idx);
    free(ix);
}

/* Allocate a table sized so that `hint` entries stay at load factor <= 0.5. */
static int vidx_alloc(GraphVIdx* ix, size_t hint) {
    size_t cap = 16;
    while (cap < (hint | 1) * 2) cap <<= 1;
    ix->key = calloc(cap, sizeof(const Expr*));
    ix->idx = calloc(cap, sizeof(int));
    if (!ix->key || !ix->idx) { free(ix->key); free(ix->idx); return 0; }
    ix->mask = cap - 1;
    return 1;
}

GraphVIdx* graph_vidx_new(size_t hint) {
    GraphVIdx* ix = calloc(1, sizeof(GraphVIdx));
    if (!ix) return NULL;
    if (!vidx_alloc(ix, hint)) { free(ix); return NULL; }
    return ix;
}

/* Slot holding `v`, or the empty slot where it belongs. */
static size_t vidx_slot(const GraphVIdx* ix, const Expr* v) {
    size_t s = (size_t)expr_hash(v) & ix->mask;
    while (ix->key[s] && !expr_eq(ix->key[s], v)) s = (s + 1) & ix->mask;
    return s;
}

int graph_vidx_get(const GraphVIdx* ix, const Expr* v) {
    size_t s = vidx_slot(ix, v);
    return ix->key[s] ? ix->idx[s] : -1;
}

/* Double the table and reinsert. Keys are borrowed, so this only moves slots. */
static int vidx_grow(GraphVIdx* ix) {
    GraphVIdx bigger;
    bigger.count = ix->count;
    if (!vidx_alloc(&bigger, (ix->mask + 1) * 2)) return 0;
    for (size_t s = 0; s <= ix->mask; s++) {
        if (!ix->key[s]) continue;
        size_t d = vidx_slot(&bigger, ix->key[s]);
        bigger.key[d] = ix->key[s];
        bigger.idx[d] = ix->idx[s];
    }
    free(ix->key);
    free(ix->idx);
    ix->key = bigger.key;
    ix->idx = bigger.idx;
    ix->mask = bigger.mask;
    return 1;
}

int graph_vidx_put(GraphVIdx* ix, const Expr* v, int index) {
    size_t s = vidx_slot(ix, v);
    if (ix->key[s]) return 0;                             /* already present */
    ix->key[s] = v;
    ix->idx[s] = index;
    ix->count++;
    if (ix->count * 2 > ix->mask + 1) {
        /* Keep the load factor bounded. A failed grow leaves the table valid,
         * just fuller -- correctness does not depend on the resize. */
        (void)vidx_grow(ix);
    }
    return 1;
}

/* Index every vertex of a canonical vertex List by position. A repeated vertex
 * keeps its FIRST index, matching the linear graph_vertex_index scan. */
static GraphVIdx* vidx_build(const Expr* verts) {
    size_t n = verts->data.function.arg_count;
    GraphVIdx* ix = graph_vidx_new(n);
    if (!ix) return NULL;
    for (size_t i = 0; i < n; i++)
        graph_vidx_put(ix, verts->data.function.args[i], (int)i);
    return ix;
}

/* ---- Edge-key set -------------------------------------------------------- *
 * Parallel-edge detection in O(E). Two normalized edges are "parallel" when
 * they connect the same endpoints in a way the graph cannot distinguish:
 *   - directed:   same head and the same ordered   pair (u, v);
 *   - undirected: same head and the same unordered pair {u, v}.
 * Directed a->b and b->a are distinct; that is allowed.
 *
 * Once both endpoints are vertex indices, an edge collapses to a 64-bit key --
 * the ordered pair when directed, the sorted pair when undirected -- so the old
 * pairwise expr_eq comparison becomes a single hashed insert. Insert returns 0
 * when the key is already present, i.e. exactly when the old test found a
 * parallel edge. */
typedef struct { uint64_t k; unsigned char directed; unsigned char used; } EKSlot;
typedef struct { EKSlot* slot; size_t mask; } EKSet;

static int ekset_init(EKSet* t, size_t ne) {
    size_t cap = 16;
    while (cap < (ne | 1) * 2) cap <<= 1;
    t->slot = calloc(cap, sizeof(EKSlot));
    if (!t->slot) return 0;
    t->mask = cap - 1;
    return 1;
}

static int ekset_insert(EKSet* t, int ia, int ib, int directed) {
    uint64_t a = (uint32_t)ia, b = (uint32_t)ib;
    if (!directed && a > b) { uint64_t tmp = a; a = b; b = tmp; }
    uint64_t k = (a << 32) | b;
    /* Mix, so that consecutive vertex indices do not probe in long runs. */
    size_t s = (size_t)((k * 0x9E3779B97F4A7C15ULL) >> 32) & t->mask;
    while (t->slot[s].used) {
        if (t->slot[s].k == k && t->slot[s].directed == (unsigned char)directed)
            return 0;                                        /* parallel edge */
        s = (s + 1) & t->mask;
    }
    t->slot[s].used = 1;
    t->slot[s].k = k;
    t->slot[s].directed = (unsigned char)directed;
    return 1;
}

int graph_vertex_index(const Expr* verts, const Expr* v) {
    if (!graph_is_list(verts)) return -1;
    for (size_t i = 0; i < verts->data.function.arg_count; i++) {
        if (expr_eq(verts->data.function.args[i], v)) return (int)i;
    }
    return -1;
}

/* True iff `opt` is a well-formed EdgeWeight -> List[n] rule, where n equals
 * `edge_count`. Shape only -- does not inspect the individual weight values,
 * which may be any expression (numeric weights are the expected case, but
 * nothing here requires it, matching how vertices are already arbitrary
 * expressions). */
static int graph_edge_weight_rule_ok(const Expr* opt, size_t edge_count) {
    if (!head_is_sym(opt, SYM_Rule) || opt->data.function.arg_count != 2) return 0;
    const Expr* key = opt->data.function.args[0];
    const Expr* val = opt->data.function.args[1];
    if (!key || key->type != EXPR_SYMBOL || key->data.symbol.name != SYM_EdgeWeight)
        return 0;
    if (!graph_is_list(val)) return 0;
    return val->data.function.arg_count == edge_count;
}

/* True iff g's shape is Graph[verts, edges] (unweighted) or
 * Graph[verts, edges, EdgeWeight -> List[n]] with n == |edges| (weighted).
 * Both `graph_is_valid` and `graph_build_adj` route through this instead of
 * duplicating an `arg_count != 2` literal -- they are two independent choke
 * points (a plan-reviewer-caught defect: widening only one left the other's
 * 8 downstream builtins rejecting every weighted graph even though GraphQ
 * reported it valid), so the arity/shape check itself must be shared, not
 * just widened identically by hand in both places. Structural shape only --
 * self-loops, parallel edges, etc. are still each caller's own job. */
static int graph_shape_ok(const Expr* g) {
    if (!head_is_sym(g, SYM_Graph)) return 0;
    size_t argc = g->data.function.arg_count;
    if (argc == 2) return 1;
    if (argc != 3) return 0;
    const Expr* edges = g->data.function.args[1];
    if (!graph_is_list(edges)) return 0;
    return graph_edge_weight_rule_ok(g->data.function.args[2],
                                      edges->data.function.arg_count);
}

/* ---- Phase 5: adjacency scaffolding --------------------------------------- */

/* Validation over an already-built vertex index; defined with graph_is_valid
 * below, and shared with graph_build_adj so the index is built only once. */
static int graph_check(const Expr* g, const GraphVIdx* ix);

void graph_adj_free(GraphAdj* a) {
    if (!a) return;
    for (int i = 0; i < a->n; i++) { free(a->out[i]); free(a->in[i]); }
    free(a->out); free(a->in);
    free(a->outdeg); free(a->indeg);
    free(a);
}

GraphAdj* graph_build_adj(const Expr* g) {
    /* Validate and index in one pass: graph_is_valid would build and throw away
     * the same vertex index, and the two fill passes below need it anyway. */
    if (!graph_shape_ok(g)) return NULL;
    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    if (!graph_is_list(verts) || !graph_is_list(edges)) return NULL;

    GraphVIdx* ix = vidx_build(verts);
    if (!ix) return NULL;
    if (!graph_check(g, ix)) { graph_vidx_free(ix); return NULL; }

    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    GraphAdj* a = calloc(1, sizeof(GraphAdj));
    if (!a) { graph_vidx_free(ix); return NULL; }
    a->n = n;
    a->verts = verts;
    a->outdeg = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    a->indeg  = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    a->out    = calloc((size_t)(n > 0 ? n : 1), sizeof(int*));
    a->in     = calloc((size_t)(n > 0 ? n : 1), sizeof(int*));

    /* Pass 1: count degrees. */
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vidx_get(ix, e->data.function.args[0]);
        int ib = graph_vidx_get(ix, e->data.function.args[1]);
        a->outdeg[ia]++; a->indeg[ib]++;
        if (kind == SYM_UndirectedEdge) { a->outdeg[ib]++; a->indeg[ia]++; }
    }
    for (int i = 0; i < n; i++) {
        a->out[i] = (a->outdeg[i] > 0) ? calloc((size_t)a->outdeg[i], sizeof(int)) : NULL;
        a->in[i]  = (a->indeg[i]  > 0) ? calloc((size_t)a->indeg[i],  sizeof(int)) : NULL;
    }

    /* Pass 2: fill (reuse degree counters as write cursors). */
    int* oc = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    int* ic = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vidx_get(ix, e->data.function.args[0]);
        int ib = graph_vidx_get(ix, e->data.function.args[1]);
        a->out[ia][oc[ia]++] = ib;  a->in[ib][ic[ib]++] = ia;
        if (kind == SYM_UndirectedEdge) {
            a->out[ib][oc[ib]++] = ia;  a->in[ia][ic[ia]++] = ib;
        }
    }
    free(oc); free(ic);
    graph_vidx_free(ix);
    return a;
}

int graph_count_components(const GraphAdj* a, const char* removed, int* active_out) {
    int n = a->n;
    char* seen = calloc((size_t)(n > 0 ? n : 1), sizeof(char));
    int* stack = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    int comps = 0, active = 0;

    for (int s = 0; s < n; s++) {
        if (removed && removed[s]) continue;
        active++;
        if (seen[s]) continue;
        /* New component: DFS over underlying undirected neighbors (out + in). */
        comps++;
        int top = 0; stack[top++] = s; seen[s] = 1;
        while (top > 0) {
            int u = stack[--top];
            for (int j = 0; j < a->outdeg[u]; j++) {
                int w = a->out[u][j];
                if ((removed && removed[w]) || seen[w]) continue;
                seen[w] = 1; stack[top++] = w;
            }
            for (int j = 0; j < a->indeg[u]; j++) {
                int w = a->in[u][j];
                if ((removed && removed[w]) || seen[w]) continue;
                seen[w] = 1; stack[top++] = w;
            }
        }
    }
    free(seen); free(stack);
    if (active_out) *active_out = active;
    return comps;
}

/* The validation body, given a vertex index already built over g's vertex List.
 * Callers have checked g's outer shape (Graph head, two List arguments).
 *
 * Rejects, in the order the MVP did: an un-normalized or 3-argument edge, a
 * self-loop, an endpoint absent from the vertex list, and a parallel edge. */
static int graph_check(const Expr* g, const GraphVIdx* ix) {
    const Expr* edges = g->data.function.args[1];
    size_t ne = edges->data.function.arg_count;

    EKSet seen;
    if (!ekset_init(&seen, ne)) return 0;

    int ok = 1;
    for (size_t i = 0; i < ne && ok; i++) {
        const Expr* edge = edges->data.function.args[i];
        const char* kind = graph_edge_kind(edge);
        if (!kind) { ok = 0; break; }              /* un-normalized / 3-arg    */

        const Expr* u = edge->data.function.args[0];
        const Expr* v = edge->data.function.args[1];
        if (expr_eq(u, v)) { ok = 0; break; }      /* self-loop                */

        int ia = graph_vidx_get(ix, u);
        int ib = graph_vidx_get(ix, v);
        if (ia < 0 || ib < 0) { ok = 0; break; }   /* endpoint not a vertex    */

        if (!ekset_insert(&seen, ia, ib, kind == SYM_DirectedEdge)) ok = 0;
    }

    free(seen.slot);
    return ok;
}

int graph_is_valid(const Expr* g) {
    if (!graph_shape_ok(g))
        return 0;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    if (!graph_is_list(verts) || !graph_is_list(edges)) return 0;

    GraphVIdx* ix = vidx_build(verts);
    if (!ix) return 0;                  /* cannot index => cannot validate */
    int ok = graph_check(g, ix);
    graph_vidx_free(ix);
    return ok;
}

Expr* graph_resolve_edge_weights(const Expr* g) {
    if (!graph_is_valid(g)) return NULL;
    size_t ne = g->data.function.args[1]->data.function.arg_count;

    if (g->data.function.arg_count == 3) {
        const Expr* wlist = g->data.function.args[2]->data.function.args[1];
        Expr** ws = (ne > 0) ? calloc(ne, sizeof(Expr*)) : NULL;
        if (ne > 0 && !ws) return NULL;
        for (size_t i = 0; i < ne; i++) ws[i] = expr_copy(wlist->data.function.args[i]);
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), ws, ne);
        free(ws);
        return out;
    }

    /* Unweighted: default every edge's weight to 1, matching Wolfram Language. */
    Expr** ws = (ne > 0) ? calloc(ne, sizeof(Expr*)) : NULL;
    if (ne > 0 && !ws) return NULL;
    for (size_t i = 0; i < ne; i++) ws[i] = expr_new_integer(1);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), ws, ne);
    free(ws);
    return out;
}

/* Approximate double value of a numeric weight, for Dijkstra's internal
 * vertex-selection comparisons ONLY -- never for a returned value (see
 * graph_weight_to_double's caller: shortestpath.c reconstructs the exact
 * GraphDistance answer separately, via Plus[] over the real Expr weights).
 * A plain (not rounding-to-nearest) conversion is fine here: it only needs to
 * preserve enough precision to compare relative distances correctly, not to
 * reproduce N[expr]'s exact rounding. Returns NAN for anything not numeric --
 * callers must gate with graph_weights_usable first. */
double graph_weight_to_double(const Expr* w) {
    if (!w) return NAN;
    switch (w->type) {
        case EXPR_INTEGER: return (double)w->data.integer;
        case EXPR_REAL:    return w->data.real;
        case EXPR_BIGINT:  return mpz_get_d(w->data.bigint);
#ifdef USE_MPFR
        case EXPR_MPFR:    return mpfr_get_d(w->data.mpfr, MPFR_RNDN);
#endif
        case EXPR_FUNCTION:
            if (head_is_sym(w, SYM_Rational) && w->data.function.arg_count == 2) {
                double p = graph_weight_to_double(w->data.function.args[0]);
                double q = graph_weight_to_double(w->data.function.args[1]);
                if (!isnan(p) && !isnan(q) && q != 0.0) return p / q;
            }
            return NAN;
        default:
            return NAN;
    }
}

int graph_weights_usable(const Expr* g) {
    if (!graph_is_valid(g) || g->data.function.arg_count != 3) return 0;
    Expr* weights = graph_resolve_edge_weights(g);
    if (!weights) return 0;

    int ok = 1;
    size_t n = weights->data.function.arg_count;
    for (size_t i = 0; i < n && ok; i++) {
        const Expr* w = weights->data.function.args[i];
        /* expr_is_numeric_like also accepts Complex (numeric-component
         * Complex[re,im]); Dijkstra needs an orderable real, so reject that
         * shape explicitly rather than reusing the check unfiltered. */
        if (!expr_is_numeric_like(w) || head_is_sym(w, SYM_Complex)) { ok = 0; break; }
        double d = graph_weight_to_double(w);
        if (isnan(d) || d < 0.0) ok = 0;
    }
    expr_free(weights);
    return ok;
}

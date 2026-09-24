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
#include <string.h>
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


void graph_adj_free(GraphAdj* a) {
    if (!a) return;
    free(a->block);
    free(a->out); free(a->in);
    free(a->outdeg); free(a->indeg);
    free(a);
}

/* CSR adjacency: every out[i]/in[i] row points into ONE int block (a->block),
 * so a build is a handful of allocations rather than two per vertex, and the
 * endpoints come pre-resolved from the validated-graph memo, so the two passes
 * below are integer-only (no expr_hash/expr_eq per edge). */
GraphAdj* graph_build_adj(const Expr* g) {
    const int *eu, *ev;
    const unsigned char* edir;
    if (!graph_edge_indices(g, &eu, &ev, &edir)) return NULL;   /* validates */
    const Expr* verts = g->data.function.args[0];
    int n = (int)verts->data.function.arg_count;
    size_t ne = g->data.function.args[1]->data.function.arg_count;
    size_t nn = (size_t)(n > 0 ? n : 1);

    GraphAdj* a = calloc(1, sizeof(GraphAdj));
    if (!a) return NULL;
    a->n = n;
    a->verts = verts;
    a->outdeg = calloc(nn, sizeof(int));
    a->indeg  = calloc(nn, sizeof(int));
    a->out    = malloc(nn * sizeof(int*));
    a->in     = malloc(nn * sizeof(int*));
    if (!a->outdeg || !a->indeg || !a->out || !a->in) { graph_adj_free(a); return NULL; }

    /* Pass 1: count degrees. */
    size_t slots = 0;
    for (size_t k = 0; k < ne; k++) {
        a->outdeg[eu[k]]++; a->indeg[ev[k]]++;
        if (!edir[k]) { a->outdeg[ev[k]]++; a->indeg[eu[k]]++; }
        slots += edir[k] ? 2 : 4;
    }
    a->block = malloc((slots > 0 ? slots : 1) * sizeof(int));
    if (!a->block) { graph_adj_free(a); return NULL; }
    int* p = a->block;
    for (int i = 0; i < n; i++) { a->out[i] = p; p += a->outdeg[i]; }
    for (int i = 0; i < n; i++) { a->in[i]  = p; p += a->indeg[i]; }

    /* Pass 2: fill, in edge order (so neighbour order matches the old
     * per-vertex build exactly). */
    int* oc = calloc(nn, sizeof(int));
    int* ic = calloc(nn, sizeof(int));
    if (!oc || !ic) { free(oc); free(ic); graph_adj_free(a); return NULL; }
    for (size_t k = 0; k < ne; k++) {
        int ia = eu[k], ib = ev[k];
        a->out[ia][oc[ia]++] = ib;  a->in[ib][ic[ib]++] = ia;
        if (!edir[k]) { a->out[ib][oc[ib]++] = ia;  a->in[ia][ic[ia]++] = ib; }
    }
    free(oc); free(ic);
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
 * self-loop, an endpoint absent from the vertex list, and a parallel edge.
 *
 * On success, when `keep` is non-NULL the edge-key set built along the way is
 * handed to the caller (who frees keep->slot) instead of being discarded, and
 * the number of directed edges is written to *ndir. When `eu` is non-NULL,
 * edge i's endpoint indices and direction go to eu/ev/edir[i]. The
 * validated-graph memo below keeps all of it, so EdgeQ and the direction
 * predicates answer in O(1) and adjacency building needs no hash lookups. */
static int graph_check_keep(const Expr* g, const GraphVIdx* ix, EKSet* keep,
                            size_t* ndir, int* eu, int* ev, unsigned char* edir) {
    const Expr* edges = g->data.function.args[1];
    size_t ne = edges->data.function.arg_count;

    EKSet seen;
    if (!ekset_init(&seen, ne)) return 0;

    int ok = 1;
    size_t directed = 0;
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

        if (kind == SYM_DirectedEdge) directed++;
        if (eu) { eu[i] = ia; ev[i] = ib; edir[i] = (unsigned char)(kind == SYM_DirectedEdge); }
        if (!ekset_insert(&seen, ia, ib, kind == SYM_DirectedEdge)) ok = 0;
    }

    if (ok && keep) { *keep = seen; *ndir = directed; }
    else free(seen.slot);
    return ok;
}

/* True iff the set holds the key ekset_insert would build for (ia, ib). */
static int ekset_contains(const EKSet* t, int ia, int ib, int directed) {
    uint64_t a = (uint32_t)ia, b = (uint32_t)ib;
    if (!directed && a > b) { uint64_t tmp = a; a = b; b = tmp; }
    uint64_t k = (a << 32) | b;
    size_t s = (size_t)((k * 0x9E3779B97F4A7C15ULL) >> 32) & t->mask;
    while (t->slot[s].used) {
        if (t->slot[s].k == k && t->slot[s].directed == (unsigned char)directed)
            return 1;
        s = (s + 1) & t->mask;
    }
    return 0;
}

/* ---- Validated-graph memo --------------------------------------------------
 * Graphs are plain Exprs, so every accessor used to re-validate its argument
 * from scratch -- building a vertex hash index and an edge-key set, O(V + E) --
 * before answering even an O(1) question: VertexQ/EdgeQ/EdgeCount on a
 * 100000-vertex graph each spent ~2-6 ms there, against ~0 ms for networkx and
 * for Mathematica's atomic Graph object.
 *
 * The memo remembers the last few graphs that passed validation, together with
 * the vertex index and edge-key set validation built anyway. It is keyed on
 * the node POINTER, which is sound because each slot holds a reference
 * (expr_copy): a node with a live reference cannot be freed, so its address
 * cannot be recycled for a different graph, and a node with refcount > 1 is
 * immutable (mutators expr_unshare first -- see struct Expr). A structurally
 * equal graph at another address simply misses and is validated afresh.
 *
 * Bounded at GRAPH_MEMO_SLOTS entries, evicted round-robin, so at most that
 * many graphs are kept alive past their last user reference. Module-static
 * with no locking, following the g_qqbar_cache precedent (flint_qqbar.c).
 * Only valid graphs are memoized: an invalid one is re-checked every call. */
#define GRAPH_MEMO_SLOTS 8

typedef struct {
    Expr*      g;       /* owned reference; NULL = empty slot                 */
    GraphVIdx* ix;      /* vertex -> index; keys borrowed from g's vertex List */
    EKSet      ek;      /* every edge's key, as built by graph_check          */
    size_t     ndir;    /* number of DirectedEdges                            */
    int*       eu;      /* eu[k], ev[k]: endpoint indices of edge k           */
    int*       ev;
    unsigned char* edir;/* edir[k]: 1 if edge k is a DirectedEdge             */
    signed char prop[GRAPH_PROP_COUNT];  /* cached answers; -1 = not computed */
    Expr*      cached[GRAPH_CACHED_COUNT]; /* owned cached results, or NULL  */
} GraphMemo;

static void graph_memo_clear(GraphMemo* m) {
    graph_vidx_free(m->ix); free(m->ek.slot);
    free(m->eu); free(m->ev); free(m->edir);
    for (int i = 0; i < GRAPH_CACHED_COUNT; i++) if (m->cached[i]) expr_free(m->cached[i]);
    expr_free(m->g);
    memset(m, 0, sizeof(*m));
}

static GraphMemo g_graph_memo[GRAPH_MEMO_SLOTS];
static int g_graph_memo_next = 0;

/* Store a validated graph's artefacts in the next slot (round-robin), taking
 * ownership of ix/ek/eu/ev/edir. */
static GraphMemo* graph_memo_insert(const Expr* g, GraphVIdx* ix, EKSet ek, size_t ndir,
                                    int* eu, int* ev, unsigned char* edir) {
    GraphMemo* m = &g_graph_memo[g_graph_memo_next];
    g_graph_memo_next = (g_graph_memo_next + 1) % GRAPH_MEMO_SLOTS;
    if (m->g) graph_memo_clear(m);
    /* expr_copy only bumps the refcount; it never writes the node's structure,
     * so casting away const here is safe. */
    m->g = expr_copy((Expr*)g);
    m->ix = ix;
    m->ek = ek;
    m->ndir = ndir;
    m->eu = eu; m->ev = ev; m->edir = edir;
    for (int i = 0; i < GRAPH_PROP_COUNT; i++) m->prop[i] = -1;
    return m;
}

/* The memo entry for g, validating (and memoizing) it on a miss. NULL iff g is
 * not a valid graph (or validation could not allocate). */
/* True iff a and b are Graph nodes with the same arity whose arguments are the
 * very same nodes (pointer-equal) -- i.e. b is a fresh wrapper around a's
 * parts, which is what the evaluator hands a builtin when it re-evaluates a
 * stored graph (or the constructor's result) without changing any argument. */
static int same_graph_parts(const Expr* a, const Expr* b) {
    if (!a || !b || a->type != EXPR_FUNCTION || b->type != EXPR_FUNCTION) return 0;
    size_t n = a->data.function.arg_count;
    if (n != b->data.function.arg_count
        || a->data.function.head->type != EXPR_SYMBOL
        || b->data.function.head->type != EXPR_SYMBOL
        || a->data.function.head->data.symbol.name != b->data.function.head->data.symbol.name)
        return 0;
    for (size_t i = 0; i < n; i++)
        if (a->data.function.args[i] != b->data.function.args[i]) return 0;
    return 1;
}

static const GraphMemo* graph_memo(const Expr* g) {
    for (int i = 0; i < GRAPH_MEMO_SLOTS; i++)
        if (g_graph_memo[i].g == g) return &g_graph_memo[i];

    /* A fresh wrapper around a memoized graph's own argument nodes: re-key the
     * entry to the new node instead of re-validating (which cost ~58 ms for a
     * 5x10^5-edge graph on every use of a variable holding it, and burned a
     * second slot per constructed graph). Sound: the vertex/edge Lists are the
     * same nodes, now kept alive by the new wrapper's reference, so the index
     * keys and every cached answer still describe it exactly. */
    for (int i = 0; i < GRAPH_MEMO_SLOTS; i++) {
        GraphMemo* m = &g_graph_memo[i];
        if (m->g && same_graph_parts(m->g, g)) {
            Expr* old = m->g;
            m->g = expr_copy((Expr*)g);   /* refcount bump only; see graph_memo_insert */
            expr_free(old);
            return m;
        }
    }

    if (!graph_shape_ok(g)) return NULL;
    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    if (!graph_is_list(verts) || !graph_is_list(edges)) return NULL;
    size_t ne = edges->data.function.arg_count, nb = ne > 0 ? ne : 1;
    /* Un-normalized edges (u -> v sugar, as Graph's constructor sees on every
     * fresh call) fail here in O(E) pointer checks, before any hashing. */
    for (size_t k = 0; k < ne; k++)
        if (!graph_edge_kind(edges->data.function.args[k])) return NULL;

    GraphVIdx* ix = vidx_build(verts);
    int* eu = malloc(nb * sizeof(int));
    int* ev = malloc(nb * sizeof(int));
    unsigned char* edir = malloc(nb);
    EKSet ek;
    size_t ndir = 0;
    if (!ix || !eu || !ev || !edir
        || !graph_check_keep(g, ix, &ek, &ndir, eu, ev, edir)) {
        graph_vidx_free(ix); free(eu); free(ev); free(edir);
        return NULL;
    }
    return graph_memo_insert(g, ix, ek, ndir, eu, ev, edir);
}

int graph_memo_seed(const Expr* g, GraphVIdx* ix, int* eu, int* ev, unsigned char* edir) {
    size_t ne = g->data.function.args[1]->data.function.arg_count;
    EKSet ek;
    int ok = ekset_init(&ek, ne);
    size_t ndir = 0;
    for (size_t k = 0; k < ne && ok; k++) {
        if (eu[k] == ev[k]) ok = 0;                    /* self-loop     */
        else if (!ekset_insert(&ek, eu[k], ev[k], edir[k])) ok = 0;  /* parallel */
        ndir += edir[k];
    }
    if (!ok) {
        free(ek.slot); graph_vidx_free(ix); free(eu); free(ev); free(edir);
        return 0;
    }
    graph_memo_insert(g, ix, ek, ndir, eu, ev, edir);
    return 1;
}

int graph_vertex_position(const Expr* g, const Expr* v) {
    const GraphMemo* m = graph_memo(g);
    if (!m) return -2;
    return graph_vidx_get(m->ix, v);
}

int graph_has_edge(const Expr* g, const Expr* u, const Expr* v, int directed) {
    const GraphMemo* m = graph_memo(g);
    if (!m) return -1;
    int ia = graph_vidx_get(m->ix, u);
    int ib = graph_vidx_get(m->ix, v);
    if (ia < 0 || ib < 0) return 0;
    return ekset_contains(&m->ek, ia, ib, directed);
}

long graph_directed_edge_count(const Expr* g) {
    const GraphMemo* m = graph_memo(g);
    return m ? (long)m->ndir : -1;
}

int graph_edge_indices(const Expr* g, const int** eu, const int** ev,
                       const unsigned char** directed) {
    const GraphMemo* m = graph_memo(g);
    if (!m) return 0;
    *eu = m->eu; *ev = m->ev; *directed = m->edir;
    return 1;
}

int graph_prop_get(const Expr* g, GraphProp p) {
    const GraphMemo* m = graph_memo(g);
    return m ? m->prop[p] : -1;
}

void graph_prop_set(const Expr* g, GraphProp p, int value) {
    GraphMemo* m = (GraphMemo*)graph_memo(g);
    if (m) m->prop[p] = (signed char)(value ? 1 : 0);
}

Expr* graph_cached_get(const Expr* g, GraphCached c) {
    GraphMemo* m = (GraphMemo*)graph_memo(g);
    return (m && m->cached[c]) ? expr_copy(m->cached[c]) : NULL;
}

void graph_cached_set(const Expr* g, GraphCached c, Expr* value) {
    GraphMemo* m = (GraphMemo*)graph_memo(g);
    if (!m) return;
    if (m->cached[c]) expr_free(m->cached[c]);
    m->cached[c] = expr_copy(value);
}

int graph_is_valid(const Expr* g) {
    if (!graph_shape_ok(g))
        return 0;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    if (!graph_is_list(verts) || !graph_is_list(edges)) return 0;

    return graph_memo(g) != NULL;
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

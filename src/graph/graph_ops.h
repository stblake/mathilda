#ifndef GRAPH_OPS_H
#define GRAPH_OPS_H

/* graph_ops.h - graph editing, transforms, set operations, structural
 * predicates and cycle/path finders (the "ops" stream of src/graph/).
 *
 * Every builtin here follows the SPEC section 4 contract used across
 * src/graph/: it borrows `res`, returns a fresh Expr* on success (the evaluator
 * frees res) or NULL to leave the call unevaluated.
 *
 * Graph results are canonical Graph[List verts, List edges(, EdgeWeight -> w)]
 * expressions that are ALREADY registered with the validated-graph memo
 * (graph_memo_seed) from the integer endpoint arrays the operation computed, so
 * the constructor's re-validation on return is an O(1) memo hit and the first
 * accessor on the result needs no hashing. Vertex and edge nodes are shared
 * with the input graphs (expr_copy) wherever they survive the edit.
 *
 * Implementation files: gops_common.c (views, result builder, edge parsing),
 * gops_edit.c, gops_transform.c, gops_setops.c, gops_preds.c, gops_cycles.c,
 * gops_init.c (registration). */

#include "expr.h"
#include "graph.h"
#include <stdint.h>

/* ---- Builtins ------------------------------------------------------------- */
/* Editing (gops_edit.c) */
Expr* builtin_vertex_add(Expr* res);          /* VertexAdd[g, v | {v..}]         */
Expr* builtin_vertex_delete(Expr* res);       /* VertexDelete[g, v | {v..}]      */
Expr* builtin_edge_add(Expr* res);            /* EdgeAdd[g, e | {e..}]           */
Expr* builtin_edge_delete(Expr* res);         /* EdgeDelete[g, e | {e..}]        */
Expr* builtin_subgraph(Expr* res);            /* Subgraph[g, vlist]              */
Expr* builtin_neighborhood_graph(Expr* res);  /* NeighborhoodGraph[g, v(, k)]    */
Expr* builtin_vertex_replace(Expr* res);      /* VertexReplace[g, rules]         */
Expr* builtin_edge_rules(Expr* res);          /* EdgeRules[g]                    */
Expr* builtin_vertex_index(Expr* res);        /* VertexIndex[g, v]               */
Expr* builtin_edge_index(Expr* res);          /* EdgeIndex[g, e]                 */
Expr* builtin_index_graph(Expr* res);         /* IndexGraph[g(, r)]              */

/* Transforms (gops_transform.c) */
Expr* builtin_graph_complement(Expr* res);    /* GraphComplement[g]              */
Expr* builtin_reverse_graph(Expr* res);       /* ReverseGraph[g]                 */
Expr* builtin_undirected_graph(Expr* res);    /* UndirectedGraph[g]              */
Expr* builtin_directed_graph(Expr* res);      /* DirectedGraph[g(, "Acyclic")]   */
Expr* builtin_line_graph(Expr* res);          /* LineGraph[g]                    */

/* Set operations (gops_setops.c) */
Expr* builtin_graph_union(Expr* res);         /* GraphUnion[g1, g2, ...]         */
Expr* builtin_graph_intersection(Expr* res);  /* GraphIntersection[g1, g2, ...]  */
Expr* builtin_graph_difference(Expr* res);    /* GraphDifference[g1, g2]         */
Expr* builtin_graph_disjoint_union(Expr* res);/* GraphDisjointUnion[g1, g2, ...] */

/* Predicates (gops_preds.c) */
Expr* builtin_simple_graph_q(Expr* res);      /* SimpleGraphQ[g]                 */
Expr* builtin_loop_free_graph_q(Expr* res);   /* LoopFreeGraphQ[g]               */
Expr* builtin_mixed_graph_q(Expr* res);       /* MixedGraphQ[g]                  */
Expr* builtin_weighted_graph_q(Expr* res);    /* WeightedGraphQ[g]               */
Expr* builtin_edge_weighted_graph_q(Expr* res);/* EdgeWeightedGraphQ[g]          */
Expr* builtin_path_graph_q(Expr* res);        /* PathGraphQ[g]                   */
Expr* builtin_eulerian_graph_q(Expr* res);    /* EulerianGraphQ[g]               */

/* Cycles and paths (gops_cycles.c) */
Expr* builtin_find_eulerian_cycle(Expr* res); /* FindEulerianCycle[g(, n)]       */
Expr* builtin_find_cycle(Expr* res);          /* FindCycle[g(, k(, n))]          */
Expr* builtin_find_path(Expr* res);           /* FindPath[g, s, t(, k(, n))]     */

/* Registers every builtin above (gops_init.c); called at the end of graph_init. */
void graph_ops_init(void);

/* ---- Shared internals (gops_common.c) -------------------------------------
 * Not builtins; used only by the gops_*.c translation units. */

/* A read-only view of a valid graph: its parts as raw arrays plus the memo's
 * integer endpoints. eu/ev/edir are BORROWED from the memo (valid until another
 * graph call may evict g) unless gops_view_own() copied them. */
typedef struct {
    const Expr*          g;
    Expr* const*         verts;   /* nv vertex nodes (borrowed)               */
    size_t               nv;
    Expr* const*         edges;   /* ne edge nodes (borrowed)                 */
    size_t               ne;
    Expr* const*         weights; /* ne weight nodes, or NULL when unweighted */
    const int*           eu;
    const int*           ev;
    const unsigned char* edir;
    size_t               ndir;    /* number of directed edges                 */
    int                  owned;   /* eu/ev/edir are private copies            */
} GopsView;

/* Fill *v for g. Returns 0 (v untouched) if g is not a valid graph. */
int  gops_view(const Expr* g, GopsView* v);
/* Replace the borrowed endpoint arrays with private copies, so the view
 * survives later graph calls (needed when several graphs are open at once).
 * Returns 0 on allocation failure. */
int  gops_view_own(GopsView* v);
void gops_view_free(GopsView* v);

/* Shared edge-head symbol nodes, created on first use; free with
 * gops_heads_free. Building one head per edge (interning per edge) is
 * measurably slow at 10^5 edges. */
typedef struct { Expr* h[2]; } GopsHeads;      /* h[0] Undirected, h[1] Directed */
/* A fresh edge node dir ? DirectedEdge[u,v] : UndirectedEdge[u,v]; takes
 * ownership of u and v. */
Expr* gops_edge(GopsHeads* hs, int directed, Expr* u, Expr* v);
void  gops_heads_free(GopsHeads* hs);

/* Classify a user-supplied edge: 1 = directed (DirectedEdge[u,v] or u->v),
 * 0 = undirected (UndirectedEdge[u,v] or u<->v), -1 = not a 2-argument edge.
 * On success *u / *v borrow the endpoints. */
int gops_parse_edge(const Expr* e, const Expr** u, const Expr** v);

/* Assemble Graph[List verts, List edges(, EdgeWeight -> List weights)] and
 * register it with the validated-graph memo. Takes ownership of every array
 * (the Expr* contents move into the result; the C arrays are freed), including
 * weights (NULL = unweighted) and eu/ev/edir. The vertex nodes must be
 * pairwise distinct (SameQ); eu/ev/edir[k] describe edges[k] by vertex
 * position. Returns NULL -- having freed everything -- if the edges contain a
 * self-loop or a parallel pair (the result would not be a valid graph) or on
 * allocation failure. */
Expr* gops_graph_new(Expr** verts, size_t nv, Expr** edges, size_t ne,
                     Expr** weights, int* eu, int* ev, unsigned char* edir);

/* Items of a "v | {v..}" argument held in *slot: a List's elements, or the
 * single argument itself (the slot). Borrowed; the count goes to *n. */
Expr* const* gops_items(Expr* const* slot, size_t* n);

/* 1 if the viewed graph has an Eulerian cycle, 0 if not, -1 if it is mixed
 * (undecided), -2 on allocation failure (gops_preds.c). */
int gops_eulerian(const GopsView* v);

/* Checked allocation helpers: calloc/malloc of n (>= 1) elements. */
void* gops_calloc(size_t n, size_t size);
void* gops_malloc(size_t n, size_t size);

/* The symbol True/False. */
Expr* gops_truth(int b);

/* Fresh List node over items (moves the Expr* contents; frees the array). */
Expr* gops_list_take(Expr** items, size_t n);

/* Incidence lists in CSR form: vertex i's entries are start[i]..start[i+1]-1,
 * each a neighbour nbr[] reached through edge eid[], in EdgeList order.
 *   GOPS_INC_ALL  every edge in both endpoints' lists (direction ignored);
 *   GOPS_INC_OUT  arcs usable forwards: directed u->v in u's list only,
 *                 undirected edges in both lists;
 *   GOPS_INC_IN   the reverse: directed u->v in v's list only. */
enum { GOPS_INC_ALL = 0, GOPS_INC_OUT = 1, GOPS_INC_IN = 2,
       GOPS_INC_OUT_NBR = 3 /* GOPS_INC_OUT without eid[] (NULL) */ };
typedef struct { int n; int* start; int* nbr; int* eid; } GopsInc;
int  gops_inc_build(const GopsView* v, int mode, GopsInc* inc);   /* 0 on OOM */
void gops_inc_free(GopsInc* inc);
/* The same CSR from a small per-graph cache (one slot per mode, holding a
 * reference to the graph). Borrowed: valid until the next call with that mode
 * for a different graph. NULL on OOM. */
const GopsInc* gops_inc_cached(const GopsView* v, int mode);

/* Stable counting sort of perm[0..m-1] (item ids) by key[item] in [0, nkeys).
 * Chain calls least-significant key first for a multi-key order. 0 on OOM. */
int gops_csort(int* perm, size_t m, const int* key, int nkeys);

/* True iff e contains a pattern construct (Blank, Pattern, PatternTest,
 * Condition, Alternatives, Except, ...), i.e. is meant as a pattern. */
int gops_has_pattern(const Expr* e);
/* MatchQ[e, patt]. */
int gops_matchq(Expr* e, Expr* patt);
/* True iff e is Infinity (the symbol or DirectedInfinity[1]). */
int gops_is_infinity(const Expr* e);

/* Open-addressed map from a 64-bit edge key to an int, sized for n keys (it
 * does not grow). gops_edge_key packs endpoint indices and direction; an
 * undirected key is orientation-free. put() inserts unless present and returns
 * the stored value either way; get() returns -1 when absent. */
typedef struct { uint64_t* key; int* val; size_t mask; } GopsKeySet;
uint64_t gops_edge_key(int a, int b, int directed);
int  gops_keyset_init(GopsKeySet* s, size_t n);             /* 0 on OOM */
void gops_keyset_free(GopsKeySet* s);
int  gops_keyset_put(GopsKeySet* s, uint64_t k, int val);
int  gops_keyset_get(const GopsKeySet* s, uint64_t k);

#endif /* GRAPH_OPS_H */

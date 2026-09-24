#ifndef GRAPH_H
#define GRAPH_H

#include "expr.h"

/* Graph subsystem entry points, registered by graph_init().
 *
 * Graphs are represented as ordinary Expr trees -- no new EXPR_* tag:
 *
 *     Graph[ List[v1, v2, ...], List[edge1, edge2, ...] ]
 *     Graph[ List[v1, v2, ...], List[edge1, edge2, ...], EdgeWeight -> List[w1, ...] ]
 *
 * where each edge is DirectedEdge[u, v] or UndirectedEdge[u, v]. Rule/->
 * and TwoWayRule/<-> are accepted as parse-time sugar and normalized on
 * construction. Vertices are arbitrary expressions. The optional 3rd argument
 * attaches a weight to each edge, matched by position (weights[i] belongs to
 * edges[i]); a graph without it is unweighted, and every accessor treats an
 * unweighted edge's weight as 1 (see graph_resolve_edge_weights).
 *
 * This mirrors the src/linalg/ layout: one builtin per translation unit,
 * with the builtin_* prototypes declared here and registered in graph.c.
 * The set of prototypes grows phase by phase as builtins land.
 *
 * Ownership follows the SPEC §4 contract: a builtin takes ownership of its
 * argument `res`, returns a new Expr* on success (the evaluator frees res),
 * or NULL to leave the expression unevaluated (the evaluator retains res).
 */

/* Graph[...] construction: normalizes edge sugar (Rule/TwoWayRule ->
 * Directed/UndirectedEdge), derives vertices when omitted, validates, and
 * returns the canonical Graph[List[verts], List[edges]]. Returns NULL to
 * leave malformed input (self-loops, parallel edges, 3-arg edges, unknown
 * edge endpoints) unevaluated, and also NULL when the argument is already
 * canonical (so evaluation reaches a fixed point). */
Expr* builtin_graph(Expr* res);

/* GraphQ[g]: True iff g is a canonical, valid graph; False otherwise. */
Expr* builtin_graph_q(Expr* res);

/* Module initializer: registers builtins, attributes, and docstrings.
 * Called from core_init() in src/core.c. */
void graph_init(void);

/* ---- Shared helpers (src/graph/graph_util.c) ------------------------------
 * Used across the subsystem (constructor, predicates, printer, and later the
 * query/matrix/algorithm builtins). All are read-only over their arguments. */

/* True iff e is a List[...] function node. */
int graph_is_list(const Expr* e);

/* If e is a normalized 2-argument edge, returns the interned edge-head pointer
 * (SYM_DirectedEdge or SYM_UndirectedEdge); otherwise NULL. */
const char* graph_edge_kind(const Expr* e);

/* True iff g is a canonical, valid graph: Graph[List verts, List edges] where
 * every edge is a 2-arg DirectedEdge/UndirectedEdge, there are no self-loops,
 * no parallel/duplicate edges, and every edge endpoint appears in verts. */
int graph_is_valid(const Expr* g);

/* Index of vertex v within List `verts` (linear expr_eq scan), or -1. Fine for a
 * single lookup; for a whole pass over the edges, build a GraphVIdx instead. */
int graph_vertex_index(const Expr* verts, const Expr* v);

/* ---- Validated-graph memo queries (src/graph/graph_util.c) ----------------
 * O(1) after the first call on a given graph node: validation results, the
 * vertex index and the edge-key set are memoized per node (see the memo
 * comment in graph_util.c for why pointer keying is sound). */

/* Position of v in g's VertexList (first occurrence, SameQ), -1 if v is not a
 * vertex, or -2 if g is not a valid graph. */
int graph_vertex_position(const Expr* g, const Expr* v);

/* 1 if g has the edge u->v (directed != 0: a DirectedEdge[u, v]) or u<->v
 * (directed == 0: an UndirectedEdge in either orientation), 0 if not, -1 if g
 * is not a valid graph. */
int graph_has_edge(const Expr* g, const Expr* u, const Expr* v, int directed);

/* Number of DirectedEdges in g, or -1 if g is not a valid graph. */
long graph_directed_edge_count(const Expr* g);

/* Borrowed views of g's edges as vertex indices: edge k runs eu[k] -> ev[k]
 * and is a DirectedEdge iff directed[k]. Valid until the next graph call that
 * could evict g from the memo -- copy what must outlive that. Returns 0 (and
 * leaves the outputs untouched) if g is not a valid graph. */
int graph_edge_indices(const Expr* g, const int** eu, const int** ev,
                       const unsigned char** directed);

/* Per-graph cache of computed answers, stored on g's memo entry, so a repeated
 * query on the same graph node is O(1) -- the same design as Mathematica's
 * atomic Graph object caching its properties. Safe for the same reason the
 * memo is: the node it is keyed on is kept alive and hence immutable. */
typedef enum {
    GRAPH_PROP_ACYCLIC, GRAPH_PROP_TREE, GRAPH_PROP_BIPARTITE,
    GRAPH_PROP_COUNT
} GraphProp;

/* The cached 0/1 answer, or -1 if not yet computed (or g not valid). */
int  graph_prop_get(const Expr* g, GraphProp p);
void graph_prop_set(const Expr* g, GraphProp p, int value);

typedef enum { GRAPH_CACHED_TOPOSORT, GRAPH_CACHED_COUNT } GraphCached;

/* A new reference to the cached result (caller frees), or NULL if none. */
Expr* graph_cached_get(const Expr* g, GraphCached c);
/* Stores a new reference to value (the caller keeps its own). */
void  graph_cached_set(const Expr* g, GraphCached c, Expr* value);

/* ---- Vertex index (src/graph/graph_util.c) --------------------------------
 * A hash index from vertex expression to an int, keyed on expr_hash/expr_eq.
 * It exists so that passes over the edge list -- validation, adjacency building,
 * deriving a vertex list from edges -- resolve endpoints in O(1) rather than by
 * a linear expr_eq scan, which is what made those passes O(E*V).
 *
 * Keys are BORROWED: the index stores the vertex pointers, so the expressions
 * must outlive it. Nothing is copied and nothing is freed but the table. */
typedef struct GraphVIdx GraphVIdx;

/* Table sized so `hint` entries fit without a resize. NULL on allocation
 * failure. Grows automatically if more are added. */
GraphVIdx* graph_vidx_new(size_t hint);
void       graph_vidx_free(GraphVIdx* ix);

/* Value stored for `v`, or -1 when absent. */
int graph_vidx_get(const GraphVIdx* ix, const Expr* v);

/* Insert `v` with value `index` if no expr_eq-equal key is present. Returns 1 if
 * inserted, 0 if it was already there (whose value is left untouched, so the
 * first insert of a repeated vertex wins). */
int graph_vidx_put(GraphVIdx* ix, const Expr* v, int index);

/* Constructor fast path: memoize g from an index and endpoint arrays its
 * builder already has, instead of re-deriving them by validation. The caller
 * guarantees g has the canonical shape with normalized edges, that every
 * endpoint resolved (eu/ev[k] >= 0), and that ix's keys are nodes g keeps alive
 * and ix maps each vertex to its FIRST VertexList position. This checks the
 * rest -- self-loops (eu == ev) and parallel edges -- and returns 1 (g is now
 * memoized) or 0 (g is invalid). Takes ownership of ix/eu/ev/edir either way. */
int graph_memo_seed(const Expr* g, GraphVIdx* ix, int* eu, int* ev, unsigned char* edir);

/* ---- Phase 2: query / representation builtins ----------------------------- */
Expr* builtin_vertex_list(Expr* res);      /* VertexList[g]                    */
Expr* builtin_edge_list(Expr* res);        /* EdgeList[g]                      */
Expr* builtin_vertex_count(Expr* res);     /* VertexCount[g]                   */
Expr* builtin_edge_count(Expr* res);       /* EdgeCount[g]                     */
Expr* builtin_adjacency_list(Expr* res);   /* AdjacencyList[g] / [g,v]         */
Expr* builtin_vertex_degree(Expr* res);    /* VertexDegree[g] / [g,v]          */
Expr* builtin_vertex_in_degree(Expr* res); /* VertexInDegree[g] / [g,v]        */
Expr* builtin_vertex_out_degree(Expr* res);/* VertexOutDegree[g] / [g,v]       */
Expr* builtin_directed_graph_q(Expr* res); /* DirectedGraphQ[g]                */
Expr* builtin_edge_weight(Expr* res);      /* EdgeWeight[g]                    */

/* ---- Phase 3: matrix views (linalg interop) ------------------------------- */
Expr* builtin_adjacency_matrix(Expr* res); /* AdjacencyMatrix[g]               */
Expr* builtin_incidence_matrix(Expr* res); /* IncidenceMatrix[g]               */
Expr* builtin_adjacency_graph(Expr* res);  /* AdjacencyGraph[m]                */
Expr* builtin_weighted_adjacency_matrix(Expr* res); /* WeightedAdjacencyMatrix[g] */

/* Resolves g's per-edge weights in EdgeList order: a copy of the EdgeWeight
 * list when g carries one (3-arg canonical form), else List[1, 1, ..., 1]
 * (one per edge). NULL if g is not a valid graph. Shared by builtin_edge_weight
 * and builtin_weighted_adjacency_matrix so the two builtins can never disagree
 * on what "unweighted" defaults to. */
Expr* graph_resolve_edge_weights(const Expr* g);

/* True iff g is a valid, weighted (3-arg EdgeWeight) graph and every one of
 * its weights is a non-negative, non-Complex number (per expr_is_numeric_like:
 * Integer, BigInt, Real, MPFR, or Rational) -- i.e. usable by a weighted
 * shortest-path algorithm. False for an unweighted graph, a symbolic weight,
 * a negative weight, or a Complex weight: shortestpath.c falls back to
 * unweighted BFS in every one of those cases rather than erroring. */
int graph_weights_usable(const Expr* g);

/* Approximate double value of a numeric weight Expr (Integer/BigInt/Real/MPFR/
 * Rational), for ranking/comparison use only -- NAN if w isn't one of those
 * shapes. Callers needing an exact returned value (e.g. GraphDistance) must
 * reconstruct it separately via real Expr arithmetic, not from this. */
double graph_weight_to_double(const Expr* w);

/* ---- Phase 4: graph generators -------------------------------------------- */
Expr* builtin_complete_graph(Expr* res);   /* CompleteGraph[n]                 */
Expr* builtin_cycle_graph(Expr* res);      /* CycleGraph[n]                    */
Expr* builtin_path_graph(Expr* res);       /* PathGraph[n] / PathGraph[{...}]  */
Expr* builtin_star_graph(Expr* res);       /* StarGraph[n]                     */
Expr* builtin_random_graph(Expr* res);     /* RandomGraph[{n,m}] / [{n,m},k]   */

/* ---- Phase 5: shared adjacency scaffolding (graph_util.c) ------------------
 * Integer-indexed adjacency derived from a validated graph. Vertex i is
 * verts[i] (canonical order). `out`/`in` hold successor/predecessor indices;
 * an UndirectedEdge{a,b} contributes symmetrically to both, so the underlying
 * undirected neighborhood of v is out[v] together with in[v]. Built on demand
 * per algorithm call (linear expr_eq indexing; documented upgrade path is an
 * expr_hash map). The caller owns the result and frees it with graph_adj_free. */
typedef struct GraphAdj {
    int   n;              /* number of vertices                                */
    const Expr* verts;    /* borrowed: the vertex List of the source graph     */
    int*  outdeg; int** out;   /* successors:   out[i][0..outdeg[i]-1]         */
    int*  indeg;  int** in;    /* predecessors: in[i][0..indeg[i]-1]           */
    int*  block;               /* CSR storage every out[i]/in[i] points into   */
} GraphAdj;

GraphAdj* graph_build_adj(const Expr* g);   /* NULL if g is not a valid graph  */
void      graph_adj_free(GraphAdj* a);

/* Count connected components of the underlying undirected graph, considering
 * only vertices with removed[i]==0 (removed may be NULL = none removed). Writes
 * the number of active vertices to *active_out when non-NULL. */
int graph_count_components(const GraphAdj* a, const char* removed, int* active_out);

/* ---- Phase 5: search & computation builtins ------------------------------- */
Expr* builtin_find_shortest_path(Expr* res); /* FindShortestPath[g,s,t]        */
Expr* builtin_graph_distance(Expr* res);     /* GraphDistance[g,s,t]           */
Expr* builtin_connected_components(Expr* res);          /* ConnectedComponents  */
Expr* builtin_weakly_connected_components(Expr* res);   /* Weakly...            */
Expr* builtin_strongly_connected_components(Expr* res); /* Strongly...          */
Expr* builtin_find_spanning_tree(Expr* res);            /* FindSpanningTree     */
Expr* builtin_connected_graph_q(Expr* res);             /* ConnectedGraphQ      */
Expr* builtin_vertex_connectivity(Expr* res);           /* VertexConnectivity   */

/* ---- Phase 5: vertex colouring internals (src/graph/vertexcoloring.c) -----
 * Non-static so they can be unit-tested directly: an exact chromatic-number
 * search has to be proven minimal before the head that answers with it is
 * registered, and a static function is not reachable from another translation
 * unit. All three are read-only over `a` and write only their out-parameters.
 *
 * `colour` buffers are caller-allocated, length a->n, and hold 1-based colour
 * numbers indexed by GraphAdj vertex index (i.e. VertexList position). */

/* DSATUR: fills colour[] with a valid -- not necessarily minimal -- colouring
 * and returns the number of colours it used, an UPPER bound on chi(g).
 * Returns 0 for the empty graph, and 0 on allocation failure. */
int fvc_dsatur_bound(const GraphAdj* a, int* colour);

/* Multi-start greedy clique: returns the size of the largest clique found, a
 * LOWER bound on chi(g). 0 for the empty graph, otherwise at least 1. */
int fvc_clique_bound(const GraphAdj* a);

/* Exact: returns chi(g) and fills colour[] with a colouring achieving it.
 * Seeds bounds from the two functions above and short-circuits when they meet,
 * otherwise refutes each smaller k by exhaustive backtracking. Writes the
 * number of backtracking nodes visited to *steps_out when non-NULL -- zero
 * exactly when the bounds met and no search was needed. */
int fvc_search(const GraphAdj* a, int* colour, long* steps_out);

/* FindVertexColoring[g]: a MINIMAL colouring, as 1-based integers in
 * VertexList order. Returns NULL (unevaluated) rather than a valid-but-larger
 * colouring whenever minimality cannot be proven -- above FVC_MAX_VERTICES, or
 * when the search exceeds FVC_MAX_STEPS. */
Expr* builtin_find_vertex_coloring(Expr* res);   /* FindVertexColoring[g]      */

/* ---- Structural predicates & ordering (graphprops.c, membership.c, acyclic.c)
 * Every *Q predicate gives False for a non-graph; TopologicalSort is left
 * unevaluated for a non-graph, a cyclic graph, or one with an undirected edge. */
Expr* builtin_undirected_graph_q(Expr* res);  /* UndirectedGraphQ[g]           */
Expr* builtin_empty_graph_q(Expr* res);       /* EmptyGraphQ[g]                */
Expr* builtin_complete_graph_q(Expr* res);    /* CompleteGraphQ[g] / [g,vlist] */
Expr* builtin_bipartite_graph_q(Expr* res);   /* BipartiteGraphQ[g]            */
Expr* builtin_vertex_q(Expr* res);            /* VertexQ[g,v]                  */
Expr* builtin_edge_q(Expr* res);              /* EdgeQ[g,e]                    */
Expr* builtin_acyclic_graph_q(Expr* res);     /* AcyclicGraphQ[g]              */
Expr* builtin_tree_graph_q(Expr* res);        /* TreeGraphQ[g]                 */
Expr* builtin_topological_sort(Expr* res);    /* TopologicalSort[g] / [{rules}]*/

/* ---- Phase 6: visualization ----------------------------------------------- */
Expr* builtin_graph_plot(Expr* res);        /* GraphPlot[g] -> Graphics[...]   */

#endif /* GRAPH_H */

#ifndef GRAPH_H
#define GRAPH_H

#include "expr.h"

/* Graph subsystem entry points, registered by graph_init().
 *
 * Graphs are represented as ordinary Expr trees -- no new EXPR_* tag:
 *
 *     Graph[ List[v1, v2, ...], List[edge1, edge2, ...] ]
 *
 * where each edge is DirectedEdge[u, v] or UndirectedEdge[u, v]. Rule/->
 * and TwoWayRule/<-> are accepted as parse-time sugar and normalized on
 * construction. Vertices are arbitrary expressions.
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

/* Index of vertex v within List `verts` (linear expr_eq scan), or -1. */
int graph_vertex_index(const Expr* verts, const Expr* v);

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

/* ---- Phase 3: matrix views (linalg interop) ------------------------------- */
Expr* builtin_adjacency_matrix(Expr* res); /* AdjacencyMatrix[g]               */
Expr* builtin_incidence_matrix(Expr* res); /* IncidenceMatrix[g]               */
Expr* builtin_adjacency_graph(Expr* res);  /* AdjacencyGraph[m]                */

/* ---- Phase 4: graph generators -------------------------------------------- */
Expr* builtin_complete_graph(Expr* res);   /* CompleteGraph[n]                 */
Expr* builtin_cycle_graph(Expr* res);      /* CycleGraph[n]                    */
Expr* builtin_path_graph(Expr* res);       /* PathGraph[n] / PathGraph[{...}]  */
Expr* builtin_random_graph(Expr* res);     /* RandomGraph[{n, m}]              */

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

/* ---- Phase 6: visualization ----------------------------------------------- */
Expr* builtin_graph_plot(Expr* res);        /* GraphPlot[g] -> Graphics[...]   */

#endif /* GRAPH_H */

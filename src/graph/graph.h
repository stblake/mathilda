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

/* Same validity check but for an arbitrary outer head (Graph or Graph3D). */
int graph_is_valid_head(const Expr* g, const char* head_sym);

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
Expr* builtin_star_graph(Expr* res);       /* StarGraph[n]                     */
Expr* builtin_wheel_graph(Expr* res);      /* WheelGraph[n]                    */
Expr* builtin_grid_graph(Expr* res);       /* GridGraph[{d1,d2,...}]           */
Expr* builtin_hypercube_graph(Expr* res);  /* HypercubeGraph[k]                */

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
Expr* builtin_bipartite_graph_q(Expr* res);             /* BipartiteGraphQ      */
Expr* builtin_vertex_eccentricity(Expr* res);           /* VertexEccentricity   */
Expr* builtin_graph_diameter(Expr* res);                /* GraphDiameter        */
Expr* builtin_graph_radius(Expr* res);                  /* GraphRadius          */
Expr* builtin_graph_center(Expr* res);                  /* GraphCenter          */
Expr* builtin_graph_periphery(Expr* res);               /* GraphPeriphery       */
Expr* builtin_acyclic_graph_q(Expr* res);               /* AcyclicGraphQ        */
Expr* builtin_topological_sort(Expr* res);              /* TopologicalSort      */
Expr* builtin_graph_complement(Expr* res);              /* GraphComplement      */
Expr* builtin_kirchhoff_matrix(Expr* res);              /* KirchhoffMatrix      */
Expr* builtin_edge_connectivity(Expr* res);             /* EdgeConnectivity     */
Expr* builtin_line_graph(Expr* res);                    /* LineGraph            */
Expr* builtin_eulerian_graph_q(Expr* res);              /* EulerianGraphQ       */
Expr* builtin_closeness_centrality(Expr* res);          /* ClosenessCentrality  */
Expr* builtin_transitive_closure(Expr* res);            /* TransitiveClosure    */
Expr* builtin_betweenness_centrality(Expr* res);        /* BetweennessCentrality*/
Expr* builtin_find_eulerian_cycle(Expr* res);           /* FindEulerianCycle    */
Expr* builtin_find_hamiltonian_cycle(Expr* res);        /* FindHamiltonianCycle */
Expr* builtin_graph_power(Expr* res);                   /* GraphPower[g, k]     */
Expr* builtin_find_cycle(Expr* res);                    /* FindCycle[g]         */
Expr* builtin_graph_distance_matrix(Expr* res);         /* GraphDistanceMatrix  */
Expr* builtin_graph_density(Expr* res);                 /* GraphDensity[g]      */
Expr* builtin_degree_centrality(Expr* res);             /* DegreeCentrality[g]  */
Expr* builtin_find_hamiltonian_path(Expr* res);         /* FindHamiltonianPath  */
Expr* builtin_kcore_components(Expr* res);              /* KCoreComponents[g,k] */
Expr* builtin_local_clustering_coefficient(Expr* res);  /* LocalClusteringCoeff */
Expr* builtin_global_clustering_coefficient(Expr* res); /* GlobalClusteringCoeff*/
Expr* builtin_mean_clustering_coefficient(Expr* res);   /* MeanClusteringCoeff  */
Expr* builtin_find_clique(Expr* res);                   /* FindClique[g]        */
Expr* builtin_find_independent_vertex_set(Expr* res);   /* FindIndependentVertexSet */
Expr* builtin_find_vertex_cover(Expr* res);             /* FindVertexCover[g]   */
Expr* builtin_graph_reciprocity(Expr* res);             /* GraphReciprocity[g]  */
Expr* builtin_chromatic_polynomial(Expr* res);          /* ChromaticPolynomial  */
Expr* builtin_chromatic_number(Expr* res);              /* ChromaticNumber[g]   */
Expr* builtin_degree_sequence(Expr* res);               /* DegreeSequence[g]    */
Expr* builtin_tree_graph_q(Expr* res);                  /* TreeGraphQ[g]        */
Expr* builtin_strongly_connected_graph_q(Expr* res);    /* StronglyConnectedGraphQ */
Expr* builtin_hamiltonian_graph_q(Expr* res);           /* HamiltonianGraphQ[g]  */
Expr* builtin_regular_graph_q(Expr* res);               /* RegularGraphQ[g]     */
Expr* builtin_complete_graph_q(Expr* res);              /* CompleteGraphQ[g]    */
Expr* builtin_graph_union(Expr* res);                   /* GraphUnion[g1, g2]   */
Expr* builtin_graph_intersection(Expr* res);            /* GraphIntersection    */
Expr* builtin_graph_difference(Expr* res);              /* GraphDifference      */
Expr* builtin_graph_reverse(Expr* res);                 /* ReverseGraph[g]      */
Expr* builtin_path_graph_q(Expr* res);                  /* PathGraphQ[g]        */
Expr* builtin_vertex_contract(Expr* res);               /* VertexContract[g,vs] */
Expr* builtin_pagerank_centrality(Expr* res);           /* PageRankCentrality   */
Expr* builtin_katz_centrality(Expr* res);               /* KatzCentrality[g,a]  */
Expr* builtin_graph_join(Expr* res);                    /* GraphJoin[g1, g2]    */
Expr* builtin_index_graph(Expr* res);                   /* IndexGraph[g] / [g,k] */
Expr* builtin_empty_graph_q(Expr* res);                 /* EmptyGraphQ[g]       */
Expr* builtin_mixed_graph_q(Expr* res);                 /* MixedGraphQ[g]       */
Expr* builtin_graph_product(Expr* res);                 /* GraphProduct[g1,g2,t] */
Expr* builtin_turan_graph(Expr* res);                   /* TuranGraph[n, r]     */
Expr* builtin_complete_kary_tree(Expr* res);            /* CompleteKaryTree[L,k] */
Expr* builtin_circulant_graph(Expr* res);               /* CirculantGraph[n,js] */
Expr* builtin_ladder_graph(Expr* res);                  /* LadderGraph[n]       */
Expr* builtin_cocktail_party_graph(Expr* res);          /* CocktailPartyGraph[n] */
Expr* builtin_kneser_graph(Expr* res);                  /* KneserGraph[n, k]    */
Expr* builtin_generalized_petersen_graph(Expr* res);    /* GeneralizedPetersenGraph */
Expr* builtin_friendship_graph(Expr* res);              /* FriendshipGraph[n]   */
Expr* builtin_vertex_coreness(Expr* res);               /* VertexCoreness[g]    */
Expr* builtin_transitive_reduction_graph(Expr* res);    /* TransitiveReductionGraph */
Expr* builtin_subgraph(Expr* res);                      /* Subgraph[g, {verts}] */
Expr* builtin_vertex_delete(Expr* res);                 /* VertexDelete[g, v]   */
Expr* builtin_edge_delete(Expr* res);                   /* EdgeDelete[g, e]     */
Expr* builtin_edge_add(Expr* res);                      /* EdgeAdd[g, e]        */
Expr* builtin_vertex_add(Expr* res);                    /* VertexAdd[g, v]      */
Expr* builtin_neighborhood_graph(Expr* res);            /* NeighborhoodGraph    */
Expr* builtin_graph_disjoint_union(Expr* res);          /* GraphDisjointUnion   */
Expr* builtin_edge_contract(Expr* res);                 /* EdgeContract[g, e]   */
Expr* builtin_find_independent_edge_set(Expr* res);     /* FindIndependentEdgeSet */
Expr* builtin_find_dominating_set(Expr* res);           /* FindDominatingSet[g] */
Expr* builtin_find_edge_cover(Expr* res);               /* FindEdgeCover[g]     */
Expr* builtin_find_vertex_coloring(Expr* res);          /* FindVertexColoring   */
Expr* builtin_graph_assortativity(Expr* res);           /* GraphAssortativity   */
Expr* builtin_incidence_list(Expr* res);                /* IncidenceList[g, v]  */
Expr* builtin_vertex_out_component(Expr* res);          /* VertexOutComponent   */
Expr* builtin_vertex_in_component(Expr* res);           /* VertexInComponent    */
Expr* builtin_antiprism_graph(Expr* res);               /* AntiprismGraph[n]    */
Expr* builtin_prism_graph(Expr* res);                   /* PrismGraph[n]        */
Expr* builtin_sunlet_graph(Expr* res);                  /* SunletGraph[n]       */
Expr* builtin_helm_graph(Expr* res);                    /* HelmGraph[n]         */
Expr* builtin_gear_graph(Expr* res);                    /* GearGraph[n]         */
Expr* builtin_edge_betweenness_centrality(Expr* res);   /* EdgeBetweennessCentrality */
Expr* builtin_dodecahedral_graph(Expr* res);            /* DodecahedralGraph[]  */
Expr* builtin_icosahedral_graph(Expr* res);             /* IcosahedralGraph[]   */

/* ---- Phase 6: visualization ----------------------------------------------- */
Expr* builtin_graph_plot(Expr* res);        /* GraphPlot[g] -> Graphics[...]   */
Expr* builtin_highlight_graph(Expr* res);   /* HighlightGraph[g, parts]        */

/* Vertex-coordinate layout (src/graph/layout.c). Fills caller-allocated x[]/y[]
 * (length = VertexCount[g]) with 2D coordinates under the named layout,
 * normalized to the [-1,1] box. `layout` is a Wolfram GraphLayout name (e.g.
 * "SpringElectricalEmbedding"); NULL / unknown falls back to circular. Fully
 * deterministic so notebooks reproduce. Read-only over g. */
void graph_compute_layout(const Expr* g, const char* layout, double* x, double* y);

/* Styling passed to graph_render(). Colors are borrowed RGBColor expressions
 * (or NULL for the built-in defaults). Highlight masks, when non-NULL, are
 * arrays as long as the vertex / edge list: a set entry draws that element in
 * the accent color and dims the rest. */
typedef struct GraphStyle {
    const char* layout;       /* layout name, or NULL for circular            */
    const Expr* vertex_color; /* RGBColor for all vertices, or NULL           */
    const Expr* edge_color;   /* RGBColor for all edges, or NULL              */
    double      vertex_size;  /* Disk radius; <= 0 uses the default           */
    int         show_labels;  /* draw Text vertex labels (default 1)          */
    const char* hi_vert;      /* per-vertex highlight mask, or NULL           */
    const char* hi_edge;      /* per-edge highlight mask, or NULL             */
} GraphStyle;

/* Render a validated graph to a Graphics[...] expression using `st` (NULL = all
 * defaults). Circular layout, blue vertices, labelled. Caller owns the result;
 * read-only over g and st. */
Expr* graph_render(const Expr* g, const GraphStyle* st);

/* Convenience wrapper: graph_render with all defaults. Used by the REPL to
 * auto-display a bare Graph result as a diagram. */
Expr* graph_default_graphics(const Expr* g);

/* ---- Graph3D: three-dimensional graphs ------------------------------------ */

/* Graph3D[...] — same construction/normalization as Graph, but the canonical
 * value has head Graph3D and auto-displays as a 3D node-link diagram. */
Expr* builtin_graph3d(Expr* res);

/* 3D vertex layout (src/graph/layout.c): fills caller-allocated x[]/y[]/z[]
 * (length = VertexCount) with coordinates in the [-1,1] cube. `layout` selects
 * a kernel (currently a 3D Fruchterman-Reingold spring, seeded on a sphere;
 * "SphericalEmbedding" places vertices on the sphere). Deterministic. */
void graph_compute_layout3d(const Expr* g, const char* layout,
                            double* x, double* y, double* z);

/* Render a valid Graph3D to a Graphics3D[...] expression (edges as 3D Lines,
 * vertices as a Point set). Caller owns the result. */
Expr* graph_render3d(const Expr* g, const GraphStyle* st);
Expr* graph_default_graphics3d(const Expr* g);

#endif /* GRAPH_H */

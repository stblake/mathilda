#ifndef GRAPH_HYPER_H
#define GRAPH_HYPER_H

#include "expr.h"
#include "graph.h"

/* graph_hyper.h - the native hypergraph subsystem (src/graph/hyp_*.c).
 *
 * REPRESENTATION. A hypergraph is a plain Expr, like Graph:
 *
 *     Hypergraph[ List[v1, ..., vn], List[e1, ..., em] ]
 *
 * where every hyperedge e_j is itself a List of vertices. Vertices are
 * arbitrary expressions and are pairwise distinct (SameQ). Hyperedges may be
 * repeated (a multi-hypergraph), may overlap or nest arbitrarily, may be empty,
 * and may repeat a vertex ({1, 1, 2}, as Wolfram-model states do).
 *
 * ORDERED VS. UNORDERED -- the documented choice. A hyperedge List is stored
 * exactly as given, so its order is preserved and available to ORDER-SENSITIVE
 * heads: HypergraphToGraph (ordered hyperedge {a, b, c} -> a->b, a->c, b->c, as
 * the Function Repository's HypergraphToGraph), EdgeList/InputForm round-trip,
 * HyperedgeSizes (arity = Length, counting repeats), and HypergraphEdgeDelete
 * (SameQ on the List). Every SET-THEORETIC head -- degrees, incidence,
 * duals, expansions, line graphs, components, distances, transversals --
 * reads a hyperedge as the SET of its distinct vertices. So a Wolfram-model
 * state {{1,2,3},{3,4}} is a valid Hypergraph edge list as-is, and the same
 * object answers both kinds of question.
 *
 * VALIDATED-HYPERGRAPH MEMO (hyp_util.c). Modelled on graph_util.c's: the last
 * HYP_MEMO_SLOTS valid hypergraphs are remembered by node pointer, each slot
 * holding a reference (expr_copy) so the node cannot be freed/recycled and is
 * immutable. A slot keeps the vertex index, the hyperedges as vertex-index CSR
 * (raw and de-duplicated), and -- built lazily on first need -- the
 * vertex->hyperedge incidence CSR. So every accessor after the first is O(1)
 * to validate, and the O(sum |e|) index build is paid once per object.
 */

/* True iff h is a canonical, valid Hypergraph (memoized; O(1) on a hit). */
int hypergraph_is_valid(const Expr* h);

/* Borrowed integer view of a valid hypergraph. Vertex i is verts[i]; hyperedge
 * j is ev[eoff[j] .. eoff[j+1]-1] (as given, repeats kept) and, as a set,
 * sv[soff[j] .. soff[j+1]-1] (distinct, first-occurrence order). When
 * incidence was requested, vertex i lies in the hyperedges
 * ve[voff[i] .. voff[i+1]-1], ascending, each once. All pointers are owned by
 * the memo and valid until the next hypergraph call that could evict h. */
typedef struct HypView {
    int n, m;
    const Expr* verts;
    const Expr* edges;
    const int *eoff, *ev;
    const int *soff, *sv;
    const int *voff, *ve;       /* NULL unless want_incidence */
    const GraphVIdx* ix;        /* vertex -> index                  */
} HypView;

/* Fills *v and returns 1 if h is valid (building the incidence CSR when
 * want_incidence), else 0. */
int hyp_view(const Expr* h, HypView* v, int want_incidence);

/* Accessor dispatch used by the shared graph accessors (vertexlist.c, ...)
 * when their argument is a Hypergraph rather than a Graph. Each takes the
 * builtin's `res` and returns the answer or NULL, per SPEC section 4. */
Expr* hyp_vertex_list(Expr* res);
Expr* hyp_edge_list(Expr* res);
Expr* hyp_vertex_count(Expr* res);
Expr* hyp_edge_count(Expr* res);
Expr* hyp_vertex_degree(Expr* res);
Expr* hyp_incidence_matrix(Expr* res);

/* Internal helpers shared by the hyp_*.c files. */
const char* hyp_sym_hypergraph(void);            /* interned "Hypergraph"    */
/* List of machine integers: a packed int64 buffer above the packing threshold,
 * else a plain List of Integers. NULL on allocation failure. */
Expr* hyp_int_list(const int64_t* xs, size_t n);

/* Builtins (one family per hyp_*.c file). */
Expr* builtin_hypergraph(Expr* res);                    /* Hypergraph            */
Expr* builtin_hypergraph_q(Expr* res);                  /* HypergraphQ           */
Expr* builtin_hyperedge_sizes(Expr* res);               /* HyperedgeSizes        */
Expr* builtin_hypergraph_rank(Expr* res);               /* HypergraphRank        */
Expr* builtin_hypergraph_corank(Expr* res);             /* HypergraphCorank      */
Expr* builtin_uniform_hypergraph_q(Expr* res);          /* UniformHypergraphQ    */
Expr* builtin_hypergraph_dual(Expr* res);               /* HypergraphDual        */
Expr* builtin_hypergraph_clique_expansion(Expr* res);   /* HypergraphCliqueExpansion */
Expr* builtin_hypergraph_star_expansion(Expr* res);     /* HypergraphStarExpansion   */
Expr* builtin_hypergraph_to_graph(Expr* res);           /* HypergraphToGraph     */
Expr* builtin_hypergraph_line_graph(Expr* res);         /* HypergraphLineGraph   */
Expr* builtin_hypergraph_connected_components(Expr* res); /* HypergraphConnectedComponents */
Expr* builtin_connected_hypergraph_q(Expr* res);        /* ConnectedHypergraphQ  */
Expr* builtin_hyperedge_connected_components(Expr* res);/* HyperedgeConnectedComponents */
Expr* builtin_hyperedge_distance(Expr* res);            /* HyperedgeDistance     */
Expr* builtin_hypergraph_distance(Expr* res);           /* HypergraphDistance    */
Expr* builtin_hypergraph_vertex_add(Expr* res);         /* HypergraphVertexAdd   */
Expr* builtin_hypergraph_vertex_delete(Expr* res);      /* HypergraphVertexDelete */
Expr* builtin_hypergraph_edge_add(Expr* res);           /* HypergraphEdgeAdd     */
Expr* builtin_hypergraph_edge_delete(Expr* res);        /* HypergraphEdgeDelete  */
Expr* builtin_subhypergraph(Expr* res);                 /* Subhypergraph         */
Expr* builtin_hypergraph_restriction(Expr* res);        /* HypergraphRestriction */
Expr* builtin_random_hypergraph(Expr* res);             /* RandomHypergraph      */
Expr* builtin_transversal_hypergraph(Expr* res);        /* TransversalHypergraph */
Expr* builtin_find_minimum_transversal(Expr* res);      /* FindMinimumTransversal */

/* Module initializer: registers every head above with attributes and
 * docstrings. Called once, at the end of graph_init(). */
void graph_hyper_init(void);

#endif /* GRAPH_HYPER_H */

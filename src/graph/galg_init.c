/* galg_init.c - registration for the algos stream (see graph_algos.h).
 *
 * graph_algos_init() is called once, at the end of graph_init(). Every head is
 * Protected, as in the Wolfram Language (ReadProtected does not exist in
 * Mathilda).
 */

#include "graph_algos.h"
#include "symtab.h"
#include "attr.h"

static void galg_reg(const char* name, Expr* (*fn)(Expr*), const char* doc) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(name, doc);
}

void graph_algos_init(void) {
    galg_reg("FindMaximumFlow", builtin_find_maximum_flow,
        "FindMaximumFlow[g, s, t] gives the value of a maximum flow from s to t (s and t may be lists of sources and sinks). FindMaximumFlow[g, s, t, \"prop\"] gives \"FlowValue\", \"FlowMatrix\" (dense n x n matrix of edge flows) or \"EdgeList\" (edges carrying flow, oriented along it). Capacities come from the EdgeCapacity -> {c1, ...} option (EdgeList order; default 1, EdgeWeight is ignored); VertexCapacity -> {c1, ...} caps the flow through each vertex. Undirected edges carry flow either way. Dinic's algorithm; exact for integer capacities.");
    galg_reg("FindMinimumCut", builtin_find_minimum_cut,
        "FindMinimumCut[g] gives {value, {part1, part2}}: a partition of the vertices minimizing the total weight of the edges from part1 to part2 (EdgeWeight when present, else 1). For directed graphs part1 is the source side. Nagamochi-Ibaraki for undirected graphs, max flows for directed ones.");
    galg_reg("FindEdgeCut", builtin_find_edge_cut,
        "FindEdgeCut[g] gives a minimum set of edges whose removal disconnects g (strongly, for directed graphs); FindEdgeCut[g, s, t] gives a minimum s-t edge cut. Uses EdgeWeight as capacity when present. Edges are returned in EdgeList order.");
    galg_reg("FindVertexCut", builtin_find_vertex_cut,
        "FindVertexCut[g] gives a minimum set of vertices whose removal disconnects the underlying undirected graph of g (n-1 vertices for a complete graph); FindVertexCut[g, s, t] gives a minimum set separating s from t ({} when they are adjacent).");
    galg_reg("EdgeConnectivity", builtin_edge_connectivity,
        "EdgeConnectivity[g] gives the minimum total weight of edges whose removal disconnects g (strongly, for directed graphs); EdgeConnectivity[g, s, t] gives the minimum s-t edge cut weight. EdgeWeight is used when present, else each edge counts 1.");
    galg_reg("FindIndependentEdgeSet", builtin_find_independent_edge_set,
        "FindIndependentEdgeSet[g] gives a maximum independent edge set (maximum matching) of g, ignoring edge direction. Hopcroft-Karp for bipartite graphs, Edmonds' blossom algorithm otherwise.");
    galg_reg("FindEdgeCover", builtin_find_edge_cover,
        "FindEdgeCover[g] gives a minimum edge cover of g: a smallest set of edges touching every vertex. Gives {} when g has an isolated vertex (no edge cover exists).");
    galg_reg("FindVertexCover", builtin_find_vertex_cover,
        "FindVertexCover[g] gives a minimum vertex cover of g (a smallest vertex set touching every edge), exact by branch and reduce. Stays unevaluated if optimality cannot be proven within the search budget.");
    galg_reg("FindIndependentVertexSet", builtin_find_independent_vertex_set,
        "FindIndependentVertexSet[g] gives {s} with s a maximum independent vertex set of g. FindIndependentVertexSet[g, k] finds a maximal independent set of at most k vertices, [g, {k}] of exactly k, [g, {kmin, kmax}] within the range; a third argument n (or All) gives up to n such sets. Exact; unevaluated if the search budget is exhausted.");
    galg_reg("IndependentVertexSetQ", builtin_independent_vertex_set_q,
        "IndependentVertexSetQ[g, vlist] gives True if no two vertices of vlist are adjacent in g.");
    galg_reg("VertexCoverQ", builtin_vertex_cover_q,
        "VertexCoverQ[g, vlist] gives True if every edge of g has an endpoint in vlist.");
    galg_reg("IndependentEdgeSetQ", builtin_independent_edge_set_q,
        "IndependentEdgeSetQ[g, elist] gives True if elist is a set of edges of g no two of which share a vertex.");
    galg_reg("EdgeCoverQ", builtin_edge_cover_q,
        "EdgeCoverQ[g, elist] gives True if elist is a set of edges of g touching every vertex of g.");
    galg_reg("FindClique", builtin_find_clique,
        "FindClique[g] gives {c} with c a maximum clique of g. FindClique[g, k] finds a maximal clique of at most k vertices, [g, {k}] of exactly k, [g, {kmin, kmax}] within the range; a third argument n (or All) gives up to n such cliques. For directed graphs a clique needs edges both ways. Exact branch and bound with colouring bounds.");
    galg_reg("FindKClique", builtin_find_k_clique,
        "FindKClique[g, k] gives {c} with c a largest k-clique of g: a vertex set whose members are pairwise within distance k.");
    galg_reg("FindHamiltonianCycle", builtin_find_hamiltonian_cycle,
        "FindHamiltonianCycle[g] gives {c} with c a Hamiltonian cycle of g as a list of edges, or {} if there is none. FindHamiltonianCycle[g, n] / [g, All] gives up to n / all Hamiltonian cycles. Exact backtracking with pruning; unevaluated if the search budget is exhausted.");
    galg_reg("FindHamiltonianPath", builtin_find_hamiltonian_path,
        "FindHamiltonianPath[g] gives a Hamiltonian path of g as a vertex list, or {} if there is none; FindHamiltonianPath[g, s, t] one from s to t. Exact; unevaluated if the search budget is exhausted.");
    galg_reg("HamiltonianGraphQ", builtin_hamiltonian_graph_q,
        "HamiltonianGraphQ[g] gives True if g has a Hamiltonian cycle. Exact; unevaluated if the search budget is exhausted.");
    galg_reg("IsomorphicGraphQ", builtin_isomorphic_graph_q,
        "IsomorphicGraphQ[g1, g2, ...] gives True if all the graphs are isomorphic. Colour refinement with individualization-refinement search.");
    galg_reg("FindGraphIsomorphism", builtin_find_graph_isomorphism,
        "FindGraphIsomorphism[g1, g2] gives {assoc}, an isomorphism from g1 to g2 as an association of vertices, or {} if the graphs are not isomorphic. FindGraphIsomorphism[g1, g2, n] / [g1, g2, All] gives up to n / all isomorphisms.");
    galg_reg("CanonicalGraph", builtin_canonical_graph,
        "CanonicalGraph[g] gives a canonical form of g on vertices 1..n: two graphs are isomorphic iff their canonical graphs are identical.");
    galg_reg("GraphAutomorphismGroup", builtin_graph_automorphism_group,
        "GraphAutomorphismGroup[g] gives the automorphism group of g as PermutationGroup[{Cycles[...], ...}] acting on vertex positions.");
    galg_reg("PlanarGraphQ", builtin_planar_graph_q,
        "PlanarGraphQ[g] gives True if g can be drawn in the plane without edge crossings. Linear-time left-right planarity test; gives False for non-graphs.");
}

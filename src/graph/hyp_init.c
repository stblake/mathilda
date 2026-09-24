/* hyp_init.c - graph_hyper_init: registers the hypergraph heads.
 *
 * Called once at the end of graph_init() (src/graph/graph.c). Every head is
 * Protected, as Graph's are; none is Listable or holds its arguments. See
 * docs/spec/builtins/hypergraphs.md for the full semantics and examples.
 */

#include "graph_hyper.h"
#include "symtab.h"
#include "attr.h"

static void reg(const char* name, BuiltinFunc f, const char* doc) {
    symtab_add_builtin(name, f);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(name, doc);
}

void graph_hyper_init(void) {
    reg("Hypergraph", builtin_hypergraph,
        "Hypergraph[{e1, e2, ...}] represents a hypergraph whose hyperedges e_i are "
        "Lists of vertices; the vertices are derived in first-appearance order. "
        "Hypergraph[{v1, ...}, {e1, ...}] gives the vertices explicitly. "
        "Hypergraph[g] converts a Graph. Hyperedges may repeat, overlap, nest, be "
        "empty, or repeat a vertex; their order is kept for ordered (Wolfram-model) "
        "use, while set-based heads read each as the set of its vertices.");
    reg("HypergraphQ", builtin_hypergraph_q,
        "HypergraphQ[h] gives True if h is a valid Hypergraph, and False otherwise.");
    reg("HyperedgeSizes", builtin_hyperedge_sizes,
        "HyperedgeSizes[h] gives the arity (Length) of each hyperedge of h, in "
        "EdgeList order.");
    reg("HypergraphRank", builtin_hypergraph_rank,
        "HypergraphRank[h] gives the largest hyperedge arity of h (0 if h has no "
        "hyperedges).");
    reg("HypergraphCorank", builtin_hypergraph_corank,
        "HypergraphCorank[h] gives the smallest hyperedge arity of h (0 if h has "
        "no hyperedges).");
    reg("UniformHypergraphQ", builtin_uniform_hypergraph_q,
        "UniformHypergraphQ[h] gives True if every hyperedge of h has the same "
        "arity. UniformHypergraphQ[h, k] gives True if h is k-uniform.");
    reg("HypergraphDual", builtin_hypergraph_dual,
        "HypergraphDual[h] gives the dual hypergraph: vertices 1..m (one per "
        "hyperedge of h) and, for each vertex v of h in VertexList order, the "
        "hyperedge of indices of the hyperedges containing v.");
    reg("HypergraphCliqueExpansion", builtin_hypergraph_clique_expansion,
        "HypergraphCliqueExpansion[h] gives the 2-section of h: the undirected "
        "Graph on VertexList[h] with u<->v whenever u and v share a hyperedge.");
    reg("HypergraphStarExpansion", builtin_hypergraph_star_expansion,
        "HypergraphStarExpansion[h] gives the incidence (star) expansion of h: the "
        "bipartite Graph on VertexList[h] and nodes Hyperedge[1], ..., Hyperedge[m], "
        "with v<->Hyperedge[j] whenever v lies in hyperedge j.");
    reg("HypergraphToGraph", builtin_hypergraph_to_graph,
        "HypergraphToGraph[h] converts h, read as an ordered hypergraph, to the "
        "directed Graph with v_a->v_b for every a<b in each hyperedge {v_1, ..., v_k} "
        "(the Function Repository's HypergraphToGraph). Self-loops are dropped and "
        "parallel edges merged. h may be a plain List of hyperedges.");
    reg("HypergraphLineGraph", builtin_hypergraph_line_graph,
        "HypergraphLineGraph[h] gives the line graph of h: vertices 1..m, i<->j when "
        "hyperedges i and j intersect. HypergraphLineGraph[h, s] gives the s-line "
        "graph, joining hyperedges that share at least s vertices.");
    reg("HypergraphConnectedComponents", builtin_hypergraph_connected_components,
        "HypergraphConnectedComponents[h] gives the connected components of h as "
        "Lists of vertices; two vertices are connected when a chain of hyperedges "
        "joins them. Components are ordered by their first vertex.");
    reg("ConnectedHypergraphQ", builtin_connected_hypergraph_q,
        "ConnectedHypergraphQ[h] gives True if h has at least one vertex and is "
        "connected. h may be a plain List of hyperedges.");
    reg("HyperedgeConnectedComponents", builtin_hyperedge_connected_components,
        "HyperedgeConnectedComponents[h] gives the components of the hyperedges of "
        "h under intersection, as Lists of 1-based hyperedge indices. "
        "HyperedgeConnectedComponents[h, s] gives the s-connected components: "
        "hyperedges with at least s vertices, joined when they share at least s.");
    reg("HyperedgeDistance", builtin_hyperedge_distance,
        "HyperedgeDistance[h, i, j] gives the length of a shortest walk of "
        "intersecting hyperedges from hyperedge i to hyperedge j, or Infinity. "
        "HyperedgeDistance[h, i, j, s] uses s-walks (consecutive hyperedges share "
        "at least s vertices). HyperedgeDistance[h, i] or [h, i, All, s] gives the "
        "distances from i to every hyperedge.");
    reg("HypergraphDistance", builtin_hypergraph_distance,
        "HypergraphDistance[h, u, v] gives the least number of hyperedges in a "
        "chain joining vertices u and v, or Infinity. HypergraphDistance[h, u] "
        "gives the distances from u to every vertex, in VertexList order.");
    reg("HypergraphVertexAdd", builtin_hypergraph_vertex_add,
        "HypergraphVertexAdd[h, v] adds the vertex v to h; HypergraphVertexAdd[h, "
        "{v1, ...}] adds several. Existing vertices are left alone.");
    reg("HypergraphVertexDelete", builtin_hypergraph_vertex_delete,
        "HypergraphVertexDelete[h, v] or [h, {v1, ...}] deletes the vertices and "
        "every hyperedge containing one of them, as VertexDelete does for a Graph.");
    reg("HypergraphEdgeAdd", builtin_hypergraph_edge_add,
        "HypergraphEdgeAdd[h, e] appends the hyperedge e (a List); "
        "HypergraphEdgeAdd[h, {e1, ...}] appends several. New vertices are added.");
    reg("HypergraphEdgeDelete", builtin_hypergraph_edge_delete,
        "HypergraphEdgeDelete[h, e] or [h, {e1, ...}] deletes every hyperedge "
        "identical (SameQ) to one given.");
    reg("Subhypergraph", builtin_subhypergraph,
        "Subhypergraph[h, {v1, ...}] gives the sub-hypergraph induced by the given "
        "vertices: those of them in h, and the hyperedges lying entirely among them.");
    reg("HypergraphRestriction", builtin_hypergraph_restriction,
        "HypergraphRestriction[h, {v1, ...}] intersects every hyperedge of h with the "
        "given vertices, dropping hyperedges that miss them (Berge's induced "
        "sub-hypergraph).");
    reg("RandomHypergraph", builtin_random_hypergraph,
        "RandomHypergraph[{n, m}, k] gives a k-uniform hypergraph on 1..n with m "
        "independent uniformly random k-subsets as hyperedges; RandomHypergraph[{n, "
        "m}, k, c] gives c of them. RandomHypergraph[{n, {e, a}}] and [{n, {{e1, a1}, "
        "...}}] draw e_i hyperedges of arity a_i with vertices chosen with "
        "replacement from 1..n. Honors SeedRandom.");
    reg("TransversalHypergraph", builtin_transversal_hypergraph,
        "TransversalHypergraph[h] gives the hypergraph of all minimal transversals "
        "(minimal hitting sets) of h, on the vertices of h. For a plain List of "
        "hyperedges the List of minimal transversals is returned. Left unevaluated "
        "if the search exceeds its budget.");
    reg("FindMinimumTransversal", builtin_find_minimum_transversal,
        "FindMinimumTransversal[h] gives a smallest set of vertices meeting every "
        "hyperedge of h (a minimum hitting set). Left unevaluated if some hyperedge "
        "is empty or the exact search exceeds its budget.");
}

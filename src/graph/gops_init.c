/* gops_init.c - registration of the graph "ops" builtins (graph_ops.h):
 * symbols, attributes and docstrings. graph_ops_init() is called once, at the
 * end of graph_init() in src/graph/graph.c. */

#include "graph_ops.h"
#include "symtab.h"
#include "attr.h"

static void reg(const char* name, BuiltinFunc f, const char* doc) {
    symtab_add_builtin(name, f);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(name, doc);
}

void graph_ops_init(void) {
    /* ---- Editing -------------------------------------------------------- */
    reg("VertexAdd", builtin_vertex_add,
        "VertexAdd[g, v] adds the vertex v to the graph g; VertexAdd[g, {v1, v2, "
        "...}] adds several. Vertices already present are ignored; new ones are "
        "appended in order. Edge weights are kept.");
    reg("VertexDelete", builtin_vertex_delete,
        "VertexDelete[g, v] removes the vertex v and its incident edges from g; "
        "VertexDelete[g, {v1, ...}] removes several (each must be a vertex of g, "
        "else the call stays unevaluated); VertexDelete[g, patt] removes every "
        "vertex matching patt. Vertex and edge order and edge weights are kept.");
    reg("EdgeAdd", builtin_edge_add,
        "EdgeAdd[g, e] adds the edge e to g; EdgeAdd[g, {e1, ...}] adds several. "
        "Edges may be written u->v / DirectedEdge[u,v] or u<->v / "
        "UndirectedEdge[u,v]; endpoints not yet in g become new vertices "
        "(appended). A new edge gets weight 1 in a weighted graph. Adding a "
        "self-loop or an edge g already has (a multigraph) is left unevaluated.");
    reg("EdgeDelete", builtin_edge_delete,
        "EdgeDelete[g, e] removes the edge e from g; EdgeDelete[g, {e1, ...}] "
        "removes several (each must be an edge of g, else the call stays "
        "unevaluated; an undirected edge matches either orientation); "
        "EdgeDelete[g, patt] removes every edge matching patt. Vertices, order "
        "and the remaining weights are kept.");
    reg("Subgraph", builtin_subgraph,
        "Subgraph[g, {v1, v2, ...}] gives the subgraph of g induced by the listed "
        "vertices (elements that are not vertices of g are ignored); "
        "Subgraph[g, patt] uses the vertices matching patt. Vertices come in the "
        "order given; edges in lower-triangular order of that vertex order. "
        "Edge weights are kept.");
    reg("NeighborhoodGraph", builtin_neighborhood_graph,
        "NeighborhoodGraph[g, v] gives the subgraph of g induced by v and its "
        "neighbours; NeighborhoodGraph[g, v, k] by the vertices within distance "
        "k of v (k a non-negative integer or Infinity), edge direction ignored. "
        "NeighborhoodGraph[g, {v1, ...}, k] uses several centres. Vertices: the "
        "centres, then each centre's new vertices in VertexList order.");
    reg("VertexReplace", builtin_vertex_replace,
        "VertexReplace[g, {v1 -> w1, ...}] replaces vertices of g according to "
        "the rules (applied to each vertex as by Replace, so patterns work). "
        "Vertices mapped together merge; a merge that would create a self-loop or "
        "a parallel edge leaves the call unevaluated. Edge weights are kept.");
    reg("EdgeRules", builtin_edge_rules,
        "EdgeRules[g] gives the edges of g as a list of rules u -> v (undirected "
        "edges too), in EdgeList order.");
    reg("VertexIndex", builtin_vertex_index,
        "VertexIndex[g, v] gives the position of the vertex v in VertexList[g]; "
        "VertexIndex[g, {v1, ...}] gives a list of positions.");
    reg("EdgeIndex", builtin_edge_index,
        "EdgeIndex[g, e] gives the position of the edge e in EdgeList[g] (an "
        "undirected edge matches either orientation; u->v and u<->v sugar are "
        "accepted); EdgeIndex[g, {e1, ...}] gives a list of positions.");
    reg("IndexGraph", builtin_index_graph,
        "IndexGraph[g] replaces the vertices of g by 1, 2, ..., n (in VertexList "
        "order); IndexGraph[g, r] by r, r+1, ..., r+n-1. Edge weights are kept.");

    /* ---- Transforms ----------------------------------------------------- */
    reg("GraphComplement", builtin_graph_complement,
        "GraphComplement[g] gives the complement of g: the same vertices, with an "
        "edge wherever g has none. For a directed or mixed g the complement is "
        "directed (u->v present iff no edge of g leads from u to v). Edges are "
        "in row-major VertexList order; weights are dropped.");
    reg("ReverseGraph", builtin_reverse_graph,
        "ReverseGraph[g] reverses every directed edge of g; undirected edges, "
        "the edge order and the weights are kept.");
    reg("UndirectedGraph", builtin_undirected_graph,
        "UndirectedGraph[g] gives the undirected graph underlying g: u->v and "
        "v->u become the single edge u<->v, whose weight is the sum of theirs. "
        "Edges are oriented and ordered by VertexList position. An undirected g "
        "is returned unchanged.");
    reg("DirectedGraph", builtin_directed_graph,
        "DirectedGraph[g] replaces each undirected edge u<->v of g by the pair "
        "u->v, v->u (weights duplicated). DirectedGraph[g, \"Acyclic\"] instead "
        "orients each undirected edge from the vertex earlier in VertexList to "
        "the later one, giving a DAG for undirected g.");
    reg("LineGraph", builtin_line_graph,
        "LineGraph[g] gives the line graph of g: one vertex per edge (numbered "
        "by EdgeList position); undirected edges sharing an endpoint are joined, "
        "and for directed g, edge i -> edge j when i ends where j starts. Mixed "
        "graphs are left unevaluated.");

    /* ---- Set operations ------------------------------------------------- */
    reg("GraphUnion", builtin_graph_union,
        "GraphUnion[g1, g2, ...] gives the graph whose vertices and edges are the "
        "unions of those of the gi (vertices in canonical order; an undirected "
        "edge equals its reversal). Weights are dropped.");
    reg("GraphIntersection", builtin_graph_intersection,
        "GraphIntersection[g1, g2, ...] gives the graph on the union of the "
        "vertex sets whose edges are those common to all the gi, in canonical "
        "order. Weights are dropped.");
    reg("GraphDifference", builtin_graph_difference,
        "GraphDifference[g1, g2] gives the graph on the union of the vertex sets "
        "whose edges are those of g1 not in g2, in canonical order. Weights are "
        "dropped.");
    reg("GraphDisjointUnion", builtin_graph_disjoint_union,
        "GraphDisjointUnion[g1, g2, ...] gives the disjoint union of the gi, "
        "with vertices relabelled 1, 2, ..., n (g1's first, in VertexList order). "
        "Weights are dropped.");

    /* ---- Predicates ----------------------------------------------------- */
    reg("SimpleGraphQ", builtin_simple_graph_q,
        "SimpleGraphQ[g] gives True if g is a graph with no self-loops or "
        "parallel edges -- every valid Mathilda graph -- and False otherwise.");
    reg("LoopFreeGraphQ", builtin_loop_free_graph_q,
        "LoopFreeGraphQ[g] gives True if g is a graph with no self-loops -- every "
        "valid Mathilda graph -- and False otherwise.");
    reg("MixedGraphQ", builtin_mixed_graph_q,
        "MixedGraphQ[g] gives True if g has both directed and undirected edges.");
    reg("WeightedGraphQ", builtin_weighted_graph_q,
        "WeightedGraphQ[g] gives True if g carries edge weights (EdgeWeight).");
    reg("EdgeWeightedGraphQ", builtin_edge_weighted_graph_q,
        "EdgeWeightedGraphQ[g] gives True if g carries edge weights (EdgeWeight).");
    reg("PathGraphQ", builtin_path_graph_q,
        "PathGraphQ[g] gives True if g is a path: connected with at least one "
        "vertex and, if undirected, every degree at most 2 and one edge fewer "
        "than vertices; if directed, every in- and out-degree at most 1. Mixed "
        "graphs are never paths.");
    reg("EulerianGraphQ", builtin_eulerian_graph_q,
        "EulerianGraphQ[g] gives True if g has a cycle using every edge exactly "
        "once: all degrees even (undirected) or in-degree = out-degree "
        "(directed), with all edges in one connected component. Left "
        "unevaluated for mixed graphs.");

    /* ---- Cycles and paths ----------------------------------------------- */
    reg("FindEulerianCycle", builtin_find_eulerian_cycle,
        "FindEulerianCycle[g] gives {c}, an Eulerian cycle c of g as a list of "
        "edges (Hierholzer's algorithm, linear time), or {} if there is none. "
        "Only FindEulerianCycle[g] and FindEulerianCycle[g, 1] are supported; "
        "mixed graphs are left unevaluated.");
    reg("FindCycle", builtin_find_cycle,
        "FindCycle[g] finds a cycle of g, as {list of edges}, or {}. FindCycle[g, "
        "k] finds a cycle of length at most k, FindCycle[g, {k}] of length "
        "exactly k, FindCycle[g, {kmin, kmax}] of length in that range; "
        "FindCycle[g, kspec, n] finds at most n (n an integer or All). Each "
        "cycle is reported once. Length-bounded searches are exhaustive and give "
        "up (unevaluated) after 5*10^7 steps; use TimeConstrained to bound them.");
    reg("FindPath", builtin_find_path,
        "FindPath[g, s, t] finds a path from s to t, as {vertex list}, or {} "
        "(depth-first, linear time). FindPath[g, s, t, k] finds a path with at "
        "most k edges, FindPath[g, s, t, {k}] exactly k, {kmin, kmax} in a range; "
        "FindPath[g, s, t, kspec, n] finds at most n paths (n or All), shortest "
        "first. Length-bounded searches give up (unevaluated) after 5*10^7 steps.");
}

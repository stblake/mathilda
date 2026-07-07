/* graph_init - the graph-subsystem module entry point.
 *
 * Registers the per-builtin symbols, attributes, and docstrings. Each builtin
 * lives in its own translation unit inside src/graph/ (construct.c, graphq.c,
 * ...), mirroring the src/linalg/ layout. graph_init() is called from
 * core_init() in src/core.c.
 *
 * Phase 1 registers Graph (construction/normalization, in construct.c) and
 * GraphQ (in graphq.c). Queries, matrix views, generators, algorithms, and
 * visualization arrive in later phases, each adding its registrations here.
 * The builtin implementations themselves live one-per-file in src/graph/.
 */

#include "graph.h"
#include "symtab.h"
#include "attr.h"

void graph_init(void) {
    /* Graph -- construction, normalization, canonicalization (construct.c). */
    symtab_add_builtin("Graph", builtin_graph);
    symtab_get_def("Graph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Graph",
        "Graph[v, e] represents a graph with vertices v and edges e. "
        "Graph[e] derives the vertices from the edge list. Edges are "
        "DirectedEdge[u,v] or UndirectedEdge[u,v]; u->v and u<->v are accepted "
        "as shorthand. Simple graphs only: no self-loops or parallel edges.");

    /* GraphQ -- validity predicate (graphq.c). */
    symtab_add_builtin("GraphQ", builtin_graph_q);
    symtab_get_def("GraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphQ",
        "GraphQ[g] gives True if g is a valid graph, and False otherwise.");

    /* ---- Phase 2: query / representation builtins ------------------------- */
    symtab_add_builtin("VertexList", builtin_vertex_list);
    symtab_get_def("VertexList")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexList",
        "VertexList[g] gives the list of vertices of the graph g.");

    symtab_add_builtin("EdgeList", builtin_edge_list);
    symtab_get_def("EdgeList")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("EdgeList",
        "EdgeList[g] gives the list of edges of the graph g.");

    symtab_add_builtin("VertexCount", builtin_vertex_count);
    symtab_get_def("VertexCount")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexCount",
        "VertexCount[g] gives the number of vertices in the graph g.");

    symtab_add_builtin("EdgeCount", builtin_edge_count);
    symtab_get_def("EdgeCount")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("EdgeCount",
        "EdgeCount[g] gives the number of edges in the graph g.");

    symtab_add_builtin("AdjacencyList", builtin_adjacency_list);
    symtab_get_def("AdjacencyList")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("AdjacencyList",
        "AdjacencyList[g] gives the adjacency list of g; AdjacencyList[g,v] "
        "gives the vertices adjacent to v (successors for directed edges).");

    symtab_add_builtin("VertexDegree", builtin_vertex_degree);
    symtab_get_def("VertexDegree")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexDegree",
        "VertexDegree[g] gives the list of vertex degrees; VertexDegree[g,v] "
        "gives the degree of vertex v.");

    symtab_add_builtin("VertexInDegree", builtin_vertex_in_degree);
    symtab_get_def("VertexInDegree")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexInDegree",
        "VertexInDegree[g] / VertexInDegree[g,v] gives in-degrees "
        "(incoming directed edges; undirected edges count for both).");

    symtab_add_builtin("VertexOutDegree", builtin_vertex_out_degree);
    symtab_get_def("VertexOutDegree")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexOutDegree",
        "VertexOutDegree[g] / VertexOutDegree[g,v] gives out-degrees "
        "(outgoing directed edges; undirected edges count for both).");

    symtab_add_builtin("DirectedGraphQ", builtin_directed_graph_q);
    symtab_get_def("DirectedGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DirectedGraphQ",
        "DirectedGraphQ[g] gives True if all edges of g are directed.");

    /* ---- Phase 3: matrix views (linalg interop) -------------------------- */
    symtab_add_builtin("AdjacencyMatrix", builtin_adjacency_matrix);
    symtab_get_def("AdjacencyMatrix")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("AdjacencyMatrix",
        "AdjacencyMatrix[g] gives the 0/1 adjacency matrix of g (symmetric for "
        "undirected graphs).");

    symtab_add_builtin("IncidenceMatrix", builtin_incidence_matrix);
    symtab_get_def("IncidenceMatrix")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("IncidenceMatrix",
        "IncidenceMatrix[g] gives the vertex-edge incidence matrix of g "
        "(oriented: -1 tail, +1 head for directed edges).");

    symtab_add_builtin("AdjacencyGraph", builtin_adjacency_graph);
    symtab_get_def("AdjacencyGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("AdjacencyGraph",
        "AdjacencyGraph[m] builds a graph on vertices 1..n from a 0/1 adjacency "
        "matrix m (undirected if m is symmetric, else directed).");

    /* ---- Phase 4: graph generators --------------------------------------- */
    symtab_add_builtin("CompleteGraph", builtin_complete_graph);
    symtab_get_def("CompleteGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CompleteGraph",
        "CompleteGraph[n] gives the complete graph K_n on n vertices.");

    symtab_add_builtin("CycleGraph", builtin_cycle_graph);
    symtab_get_def("CycleGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CycleGraph",
        "CycleGraph[n] gives the cycle graph on n vertices.");

    symtab_add_builtin("PathGraph", builtin_path_graph);
    symtab_get_def("PathGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PathGraph",
        "PathGraph[n] gives the path on n vertices; PathGraph[{v1,...}] the "
        "path over the given vertices.");

    symtab_add_builtin("RandomGraph", builtin_random_graph);
    symtab_get_def("RandomGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("RandomGraph",
        "RandomGraph[{n, m}] gives a random undirected graph with n vertices "
        "and m edges.");

    /* ---- Phase 5: search & computation algorithms ------------------------ */
    symtab_add_builtin("FindShortestPath", builtin_find_shortest_path);
    symtab_get_def("FindShortestPath")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindShortestPath",
        "FindShortestPath[g,s,t] gives a shortest path from s to t as a list of "
        "vertices ({} if none).");

    symtab_add_builtin("GraphDistance", builtin_graph_distance);
    symtab_get_def("GraphDistance")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphDistance",
        "GraphDistance[g,s,t] gives the length of a shortest path from s to t "
        "(Infinity if unreachable).");

    symtab_add_builtin("ConnectedComponents", builtin_connected_components);
    symtab_get_def("ConnectedComponents")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ConnectedComponents",
        "ConnectedComponents[g] gives the connected components of g (weak, on "
        "the underlying undirected graph).");

    symtab_add_builtin("WeaklyConnectedComponents", builtin_weakly_connected_components);
    symtab_get_def("WeaklyConnectedComponents")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("WeaklyConnectedComponents",
        "WeaklyConnectedComponents[g] gives the weakly connected components of g.");

    symtab_add_builtin("StronglyConnectedComponents", builtin_strongly_connected_components);
    symtab_get_def("StronglyConnectedComponents")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("StronglyConnectedComponents",
        "StronglyConnectedComponents[g] gives the strongly connected components "
        "of g (following edge directions).");

    symtab_add_builtin("FindSpanningTree", builtin_find_spanning_tree);
    symtab_get_def("FindSpanningTree")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindSpanningTree",
        "FindSpanningTree[g] gives a spanning tree (forest) of g as a graph.");

    symtab_add_builtin("ConnectedGraphQ", builtin_connected_graph_q);
    symtab_get_def("ConnectedGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ConnectedGraphQ",
        "ConnectedGraphQ[g] gives True if g is connected.");

    symtab_add_builtin("VertexConnectivity", builtin_vertex_connectivity);
    symtab_get_def("VertexConnectivity")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexConnectivity",
        "VertexConnectivity[g] gives the minimum number of vertices whose "
        "removal disconnects g.");

    /* ---- Phase 6: visualization ------------------------------------------ */
    symtab_add_builtin("GraphPlot", builtin_graph_plot);
    symtab_get_def("GraphPlot")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphPlot",
        "GraphPlot[g] gives a Graphics object drawing the graph g with a "
        "circular vertex layout.");
}

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

    symtab_add_builtin("StarGraph", builtin_star_graph);
    symtab_get_def("StarGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("StarGraph",
        "StarGraph[n] gives the star K_{1,n-1}: a central vertex joined to n-1 "
        "leaves.");

    symtab_add_builtin("WheelGraph", builtin_wheel_graph);
    symtab_get_def("WheelGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("WheelGraph",
        "WheelGraph[n] gives the wheel on n vertices: a cycle of n-1 rim "
        "vertices plus a hub joined to all of them (n >= 4).");

    symtab_add_builtin("GridGraph", builtin_grid_graph);
    symtab_get_def("GridGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GridGraph",
        "GridGraph[{d1,d2,...}] gives the k-dimensional grid graph on a "
        "d1 x d2 x ... lattice (cells adjacent when they differ by 1 in one "
        "coordinate).");

    symtab_add_builtin("HypercubeGraph", builtin_hypercube_graph);
    symtab_get_def("HypercubeGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("HypercubeGraph",
        "HypercubeGraph[k] gives the k-dimensional hypercube Q_k: 2^k vertices "
        "adjacent when they differ in one bit (k-regular, bipartite).");

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

    symtab_add_builtin("BipartiteGraphQ", builtin_bipartite_graph_q);
    symtab_get_def("BipartiteGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("BipartiteGraphQ",
        "BipartiteGraphQ[g] gives True if the underlying undirected graph is "
        "2-colorable (has no odd cycle), and False otherwise.");

    symtab_add_builtin("VertexEccentricity", builtin_vertex_eccentricity);
    symtab_get_def("VertexEccentricity")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexEccentricity",
        "VertexEccentricity[g,v] gives the greatest shortest-path distance from "
        "v to any vertex; VertexEccentricity[g] gives the list for all vertices "
        "(Infinity if some vertex is unreachable).");

    symtab_add_builtin("GraphDiameter", builtin_graph_diameter);
    symtab_get_def("GraphDiameter")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphDiameter",
        "GraphDiameter[g] gives the maximum vertex eccentricity (Infinity if g "
        "is not strongly connected).");

    symtab_add_builtin("GraphRadius", builtin_graph_radius);
    symtab_get_def("GraphRadius")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphRadius",
        "GraphRadius[g] gives the minimum vertex eccentricity (Infinity if no "
        "vertex reaches all others).");

    symtab_add_builtin("GraphCenter", builtin_graph_center);
    symtab_get_def("GraphCenter")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphCenter",
        "GraphCenter[g] gives the vertices whose eccentricity equals the graph "
        "radius.");

    symtab_add_builtin("GraphPeriphery", builtin_graph_periphery);
    symtab_get_def("GraphPeriphery")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphPeriphery",
        "GraphPeriphery[g] gives the vertices whose eccentricity equals the graph "
        "diameter.");

    symtab_add_builtin("AcyclicGraphQ", builtin_acyclic_graph_q);
    symtab_get_def("AcyclicGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("AcyclicGraphQ",
        "AcyclicGraphQ[g] gives True if g has no cycle: a DAG for a directed "
        "graph, a forest for an undirected one.");

    symtab_add_builtin("TopologicalSort", builtin_topological_sort);
    symtab_get_def("TopologicalSort")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("TopologicalSort",
        "TopologicalSort[g] gives a vertex ordering in which every edge points "
        "forward, or $Failed if g is not a directed acyclic graph.");

    symtab_add_builtin("GraphComplement", builtin_graph_complement);
    symtab_get_def("GraphComplement")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphComplement",
        "GraphComplement[g] gives the graph on the same vertices whose edges are "
        "exactly the non-edges of g (directed graphs stay directed).");

    symtab_add_builtin("KirchhoffMatrix", builtin_kirchhoff_matrix);
    symtab_get_def("KirchhoffMatrix")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KirchhoffMatrix",
        "KirchhoffMatrix[g] gives the graph Laplacian D - A (degree diagonal "
        "minus adjacency matrix); each row sums to 0.");

    symtab_add_builtin("EdgeConnectivity", builtin_edge_connectivity);
    symtab_get_def("EdgeConnectivity")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("EdgeConnectivity",
        "EdgeConnectivity[g] gives the minimum number of edges whose removal "
        "disconnects g (0 if g is already disconnected).");

    symtab_add_builtin("LineGraph", builtin_line_graph);
    symtab_get_def("LineGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LineGraph",
        "LineGraph[g] gives the line graph of g: its vertices are the edges of "
        "g, adjacent when they share an endpoint (head-to-tail if directed).");

    symtab_add_builtin("EulerianGraphQ", builtin_eulerian_graph_q);
    symtab_get_def("EulerianGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("EulerianGraphQ",
        "EulerianGraphQ[g] gives True if g has an Eulerian cycle: connected with "
        "all even degrees (undirected), or in-degree = out-degree everywhere "
        "(directed).");

    symtab_add_builtin("ClosenessCentrality", builtin_closeness_centrality);
    symtab_get_def("ClosenessCentrality")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ClosenessCentrality",
        "ClosenessCentrality[g] gives the list of closeness centralities "
        "c_i = (r_i-1)^2/((n-1) S_i), where r_i vertices are reachable from i at "
        "total distance S_i; larger means more central.");

    symtab_add_builtin("TransitiveClosure", builtin_transitive_closure);
    symtab_get_def("TransitiveClosure")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("TransitiveClosure",
        "TransitiveClosure[g] adds an edge u->v whenever v is reachable from u "
        "(directed); for an undirected graph each connected component becomes a "
        "complete graph.");

    symtab_add_builtin("BetweennessCentrality", builtin_betweenness_centrality);
    symtab_get_def("BetweennessCentrality")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("BetweennessCentrality",
        "BetweennessCentrality[g] gives, for each vertex, the number of shortest "
        "paths passing through it (fractional when paths tie); undirected pairs "
        "counted once.");

    symtab_add_builtin("FindEulerianCycle", builtin_find_eulerian_cycle);
    symtab_get_def("FindEulerianCycle")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindEulerianCycle",
        "FindEulerianCycle[g] gives an Eulerian cycle (a closed walk using every "
        "edge once) as a vertex list, or {} if g has none.");

    symtab_add_builtin("FindHamiltonianCycle", builtin_find_hamiltonian_cycle);
    symtab_get_def("FindHamiltonianCycle")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindHamiltonianCycle",
        "FindHamiltonianCycle[g] gives a Hamiltonian cycle (a closed walk "
        "visiting every vertex once) as a vertex list, or {} if g has none.");

    symtab_add_builtin("GraphPower", builtin_graph_power);
    symtab_get_def("GraphPower")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphPower",
        "GraphPower[g, k] gives the k-th power of g: the graph on the same "
        "vertices joining two vertices whenever g has a path of length at most k "
        "between them.");

    symtab_add_builtin("FindCycle", builtin_find_cycle);
    symtab_get_def("FindCycle")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindCycle",
        "FindCycle[g] gives a cycle in g as a list containing one list of its "
        "edges, or {} if g is acyclic.");

    symtab_add_builtin("GraphDistanceMatrix", builtin_graph_distance_matrix);
    symtab_get_def("GraphDistanceMatrix")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphDistanceMatrix",
        "GraphDistanceMatrix[g] gives the matrix whose (i,j) entry is the "
        "shortest-path distance from vertex i to vertex j (Infinity if "
        "unreachable).");

    symtab_add_builtin("GraphDensity", builtin_graph_density);
    symtab_get_def("GraphDensity")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphDensity",
        "GraphDensity[g] gives the fraction of possible edges present in g, an "
        "exact rational in [0, 1] (1 for a complete graph, 0 for an empty one).");

    symtab_add_builtin("DegreeCentrality", builtin_degree_centrality);
    symtab_get_def("DegreeCentrality")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DegreeCentrality",
        "DegreeCentrality[g] gives, for each vertex, the number of incident "
        "edges (degree; in-degree + out-degree for a directed graph).");

    symtab_add_builtin("FindHamiltonianPath", builtin_find_hamiltonian_path);
    symtab_get_def("FindHamiltonianPath")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindHamiltonianPath",
        "FindHamiltonianPath[g] gives a Hamiltonian path (a walk visiting every "
        "vertex once) as a vertex list, or {} if g has none.");

    symtab_add_builtin("KCoreComponents", builtin_kcore_components);
    symtab_get_def("KCoreComponents")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KCoreComponents",
        "KCoreComponents[g, k] gives the connected components of the k-core of g "
        "(the maximal subgraph in which every vertex has degree at least k), as "
        "a list of vertex lists.");

    symtab_add_builtin("LocalClusteringCoefficient", builtin_local_clustering_coefficient);
    symtab_get_def("LocalClusteringCoefficient")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LocalClusteringCoefficient",
        "LocalClusteringCoefficient[g] gives, for each vertex, the fraction of "
        "its neighbor pairs that are adjacent (0 for degree < 2).");

    symtab_add_builtin("GlobalClusteringCoefficient", builtin_global_clustering_coefficient);
    symtab_get_def("GlobalClusteringCoefficient")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GlobalClusteringCoefficient",
        "GlobalClusteringCoefficient[g] gives the graph transitivity: three "
        "times the number of triangles divided by the number of connected vertex "
        "triples (0 when there are none).");

    symtab_add_builtin("MeanClusteringCoefficient", builtin_mean_clustering_coefficient);
    symtab_get_def("MeanClusteringCoefficient")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MeanClusteringCoefficient",
        "MeanClusteringCoefficient[g] gives the average of the local clustering "
        "coefficients over all vertices.");

    symtab_add_builtin("FindClique", builtin_find_clique);
    symtab_get_def("FindClique")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindClique",
        "FindClique[g] gives a largest clique (a set of pairwise-adjacent "
        "vertices) as a list containing one vertex list.");

    symtab_add_builtin("FindIndependentVertexSet", builtin_find_independent_vertex_set);
    symtab_get_def("FindIndependentVertexSet")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindIndependentVertexSet",
        "FindIndependentVertexSet[g] gives a largest independent vertex set (a "
        "set of pairwise non-adjacent vertices) as a list containing one vertex "
        "list.");

    symtab_add_builtin("FindVertexCover", builtin_find_vertex_cover);
    symtab_get_def("FindVertexCover")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindVertexCover",
        "FindVertexCover[g] gives a minimum vertex cover (a smallest set of "
        "vertices touching every edge) as a vertex list.");

    symtab_add_builtin("GraphReciprocity", builtin_graph_reciprocity);
    symtab_get_def("GraphReciprocity")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphReciprocity",
        "GraphReciprocity[g] gives the fraction of arcs whose reverse is also "
        "present (1 for an undirected graph), as an exact rational.");

    symtab_add_builtin("ChromaticPolynomial", builtin_chromatic_polynomial);
    symtab_get_def("ChromaticPolynomial")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ChromaticPolynomial",
        "ChromaticPolynomial[g, k] gives the chromatic polynomial of g in k: the "
        "number of proper k-colorings (a polynomial for symbolic k).");

    symtab_add_builtin("ChromaticNumber", builtin_chromatic_number);
    symtab_get_def("ChromaticNumber")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ChromaticNumber",
        "ChromaticNumber[g] gives the least number of colors needed to color g "
        "so that adjacent vertices differ.");

    symtab_add_builtin("DegreeSequence", builtin_degree_sequence);
    symtab_get_def("DegreeSequence")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DegreeSequence",
        "DegreeSequence[g] gives the vertex degrees (in-degree + out-degree for "
        "a directed graph) sorted in non-increasing order.");

    symtab_add_builtin("TreeGraphQ", builtin_tree_graph_q);
    symtab_get_def("TreeGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("TreeGraphQ",
        "TreeGraphQ[g] gives True iff g is a tree: connected with no cycles "
        "(n-1 edges on n>=1 vertices).");

    symtab_add_builtin("StronglyConnectedGraphQ", builtin_strongly_connected_graph_q);
    symtab_get_def("StronglyConnectedGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("StronglyConnectedGraphQ",
        "StronglyConnectedGraphQ[g] gives True iff every vertex is reachable "
        "from every other following edge directions.");

    symtab_add_builtin("HamiltonianGraphQ", builtin_hamiltonian_graph_q);
    symtab_get_def("HamiltonianGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("HamiltonianGraphQ",
        "HamiltonianGraphQ[g] gives True iff g has a Hamiltonian cycle (a closed "
        "walk visiting every vertex once).");

    symtab_add_builtin("RegularGraphQ", builtin_regular_graph_q);
    symtab_get_def("RegularGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("RegularGraphQ",
        "RegularGraphQ[g] gives True iff every vertex has the same degree "
        "(equal in- and out-degrees for a directed graph).");

    symtab_add_builtin("CompleteGraphQ", builtin_complete_graph_q);
    symtab_get_def("CompleteGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CompleteGraphQ",
        "CompleteGraphQ[g] gives True iff every pair of distinct vertices in g "
        "is adjacent.");

    symtab_add_builtin("GraphUnion", builtin_graph_union);
    symtab_get_def("GraphUnion")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphUnion",
        "GraphUnion[g1, g2] gives the graph with the union of the vertices and "
        "the union of the edges of g1 and g2.");

    symtab_add_builtin("GraphIntersection", builtin_graph_intersection);
    symtab_get_def("GraphIntersection")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphIntersection",
        "GraphIntersection[g1, g2] gives the graph with the vertices common to "
        "g1 and g2 and the edges present in both.");

    symtab_add_builtin("GraphDifference", builtin_graph_difference);
    symtab_get_def("GraphDifference")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphDifference",
        "GraphDifference[g1, g2] gives the graph on g1's vertices with the edges "
        "of g1 that are not in g2.");

    symtab_add_builtin("ReverseGraph", builtin_graph_reverse);
    symtab_get_def("ReverseGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ReverseGraph",
        "ReverseGraph[g] reverses the direction of every edge of g (undirected "
        "edges are unchanged).");

    symtab_add_builtin("PathGraphQ", builtin_path_graph_q);
    symtab_get_def("PathGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PathGraphQ",
        "PathGraphQ[g] gives True iff g is a path graph (a tree with maximum "
        "degree at most 2).");

    symtab_add_builtin("VertexContract", builtin_vertex_contract);
    symtab_get_def("VertexContract")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexContract",
        "VertexContract[g, {v1, v2, ...}] merges the listed vertices into one, "
        "redirecting edges, dropping self-loops, and collapsing parallel edges.");

    symtab_add_builtin("PageRankCentrality", builtin_pagerank_centrality);
    symtab_get_def("PageRankCentrality")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PageRankCentrality",
        "PageRankCentrality[g] gives the PageRank of each vertex (damping 17/20) "
        "as an exact rational probability vector summing to 1.");

    symtab_add_builtin("KatzCentrality", builtin_katz_centrality);
    symtab_get_def("KatzCentrality")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KatzCentrality",
        "KatzCentrality[g, alpha] gives the Katz centrality of each vertex with "
        "attenuation alpha (base weight 1), solving (I - alpha A^T) x = 1 exactly.");

    symtab_add_builtin("GraphJoin", builtin_graph_join);
    symtab_get_def("GraphJoin")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphJoin",
        "GraphJoin[g1, g2] gives the join of g1 and g2: their disjoint union "
        "(vertices relabeled 1..n1+n2) plus every edge between the two blocks.");

    symtab_add_builtin("IndexGraph", builtin_index_graph);
    symtab_get_def("IndexGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("IndexGraph",
        "IndexGraph[g] renames the vertices of g to consecutive integers from 1 "
        "(IndexGraph[g, k] from k), remapping edges.");

    symtab_add_builtin("EmptyGraphQ", builtin_empty_graph_q);
    symtab_get_def("EmptyGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("EmptyGraphQ",
        "EmptyGraphQ[g] gives True iff g has no edges.");

    symtab_add_builtin("MixedGraphQ", builtin_mixed_graph_q);
    symtab_get_def("MixedGraphQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MixedGraphQ",
        "MixedGraphQ[g] gives True iff g has both directed and undirected edges.");

    symtab_add_builtin("GraphProduct", builtin_graph_product);
    symtab_get_def("GraphProduct")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphProduct",
        "GraphProduct[g1, g2, type] gives a product graph on V1 x V2 for type "
        "\"Cartesian\", \"Tensor\", \"Strong\", or \"Lexicographic\".");

    symtab_add_builtin("TuranGraph", builtin_turan_graph);
    symtab_get_def("TuranGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("TuranGraph",
        "TuranGraph[n, r] gives the Turan graph: the balanced complete r-partite "
        "graph on n vertices (the largest (r+1)-clique-free graph).");

    symtab_add_builtin("CompleteKaryTree", builtin_complete_kary_tree);
    symtab_get_def("CompleteKaryTree")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CompleteKaryTree",
        "CompleteKaryTree[L] / CompleteKaryTree[L, k] gives the complete k-ary "
        "tree with L levels (k=2 by default).");

    symtab_add_builtin("CirculantGraph", builtin_circulant_graph);
    symtab_get_def("CirculantGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CirculantGraph",
        "CirculantGraph[n, {j1, ...}] gives the circulant graph on n vertices "
        "joining vertex i to i +/- jk (mod n) for each jump jk.");

    symtab_add_builtin("LadderGraph", builtin_ladder_graph);
    symtab_get_def("LadderGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LadderGraph",
        "LadderGraph[n] gives the ladder graph: two n-vertex paths joined by n "
        "rungs (the product P_n x P_2).");

    symtab_add_builtin("CocktailPartyGraph", builtin_cocktail_party_graph);
    symtab_get_def("CocktailPartyGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CocktailPartyGraph",
        "CocktailPartyGraph[n] gives the cocktail-party graph K_{n x 2}: 2n "
        "vertices in n couples, each joined to all but its partner.");

    symtab_add_builtin("KneserGraph", builtin_kneser_graph);
    symtab_get_def("KneserGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KneserGraph",
        "KneserGraph[n, k] gives the Kneser graph: k-subsets of {1..n} as "
        "vertices, adjacent iff disjoint (K(5,2) is the Petersen graph).");

    symtab_add_builtin("GeneralizedPetersenGraph", builtin_generalized_petersen_graph);
    symtab_get_def("GeneralizedPetersenGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GeneralizedPetersenGraph",
        "GeneralizedPetersenGraph[n, k] gives GP(n,k): an outer n-cycle, an inner "
        "star polygon {n/k}, and spokes joining them (GP(5,2) is the Petersen graph).");

    symtab_add_builtin("FriendshipGraph", builtin_friendship_graph);
    symtab_get_def("FriendshipGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FriendshipGraph",
        "FriendshipGraph[n] gives the windmill graph of n triangles sharing one "
        "central vertex (F_1 is a triangle, F_2 the bowtie).");

    symtab_add_builtin("VertexCoreness", builtin_vertex_coreness);
    symtab_get_def("VertexCoreness")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexCoreness",
        "VertexCoreness[g] gives, for each vertex, its coreness: the largest k "
        "such that it lies in the k-core of g.");

    symtab_add_builtin("TransitiveReductionGraph", builtin_transitive_reduction_graph);
    symtab_get_def("TransitiveReductionGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("TransitiveReductionGraph",
        "TransitiveReductionGraph[g] gives the transitive reduction of a directed "
        "acyclic graph: the fewest-edge graph with the same reachability.");

    symtab_add_builtin("Subgraph", builtin_subgraph);
    symtab_get_def("Subgraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Subgraph",
        "Subgraph[g, {v1, ...}] gives the subgraph of g induced by the listed "
        "vertices (edges with both endpoints among them).");

    symtab_add_builtin("VertexDelete", builtin_vertex_delete);
    symtab_get_def("VertexDelete")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexDelete",
        "VertexDelete[g, v] (or VertexDelete[g, {v1, ...}]) removes the given "
        "vertices and all incident edges from g.");

    symtab_add_builtin("EdgeDelete", builtin_edge_delete);
    symtab_get_def("EdgeDelete")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("EdgeDelete",
        "EdgeDelete[g, e] (or EdgeDelete[g, {e1, ...}]) removes the given edges "
        "from g, keeping all vertices.");

    symtab_add_builtin("EdgeAdd", builtin_edge_add);
    symtab_get_def("EdgeAdd")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("EdgeAdd",
        "EdgeAdd[g, e] (or EdgeAdd[g, {e1, ...}]) adds the given edges to g, "
        "introducing any missing endpoints as new vertices.");

    symtab_add_builtin("VertexAdd", builtin_vertex_add);
    symtab_get_def("VertexAdd")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexAdd",
        "VertexAdd[g, v] (or VertexAdd[g, {v1, ...}]) adds the given vertices to "
        "g as isolated vertices, leaving edges unchanged.");

    symtab_add_builtin("NeighborhoodGraph", builtin_neighborhood_graph);
    symtab_get_def("NeighborhoodGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("NeighborhoodGraph",
        "NeighborhoodGraph[g, v] (or [g, v, k]) gives the subgraph induced by v "
        "and all vertices within distance k (k=1 by default).");

    symtab_add_builtin("GraphDisjointUnion", builtin_graph_disjoint_union);
    symtab_get_def("GraphDisjointUnion")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphDisjointUnion",
        "GraphDisjointUnion[g1, g2] gives the disjoint union of g1 and g2 "
        "(vertices relabeled 1..n1+n2, no edges between the two blocks).");

    symtab_add_builtin("EdgeContract", builtin_edge_contract);
    symtab_get_def("EdgeContract")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("EdgeContract",
        "EdgeContract[g, e] contracts edge e, merging its endpoints into one "
        "vertex, dropping the self-loop, and collapsing parallel edges.");

    symtab_add_builtin("FindIndependentEdgeSet", builtin_find_independent_edge_set);
    symtab_get_def("FindIndependentEdgeSet")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindIndependentEdgeSet",
        "FindIndependentEdgeSet[g] gives a maximum matching of g: a largest set "
        "of edges no two of which share a vertex.");

    symtab_add_builtin("FindDominatingSet", builtin_find_dominating_set);
    symtab_get_def("FindDominatingSet")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindDominatingSet",
        "FindDominatingSet[g] gives a minimum dominating set: a smallest set of "
        "vertices such that every vertex is in it or adjacent to it.");

    symtab_add_builtin("FindEdgeCover", builtin_find_edge_cover);
    symtab_get_def("FindEdgeCover")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindEdgeCover",
        "FindEdgeCover[g] gives a minimum edge cover: a smallest set of edges "
        "such that every vertex is incident to at least one ({} if none exists).");

    symtab_add_builtin("FindVertexColoring", builtin_find_vertex_coloring);
    symtab_get_def("FindVertexColoring")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FindVertexColoring",
        "FindVertexColoring[g] gives a proper coloring of g using the fewest "
        "colors, as a list of color indices (one per vertex).");

    symtab_add_builtin("GraphAssortativity", builtin_graph_assortativity);
    symtab_get_def("GraphAssortativity")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphAssortativity",
        "GraphAssortativity[g] gives the degree assortativity coefficient of g "
        "(exact rational in [-1,1]; Indeterminate for a regular graph).");

    symtab_add_builtin("IncidenceList", builtin_incidence_list);
    symtab_get_def("IncidenceList")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("IncidenceList",
        "IncidenceList[g, v] gives the list of edges of g incident to vertex v.");

    symtab_add_builtin("VertexOutComponent", builtin_vertex_out_component);
    symtab_get_def("VertexOutComponent")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexOutComponent",
        "VertexOutComponent[g, v] gives the vertices reachable from v (including "
        "v) following edge directions.");

    symtab_add_builtin("VertexInComponent", builtin_vertex_in_component);
    symtab_get_def("VertexInComponent")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VertexInComponent",
        "VertexInComponent[g, v] gives the vertices from which v is reachable "
        "(including v).");

    symtab_add_builtin("AntiprismGraph", builtin_antiprism_graph);
    symtab_get_def("AntiprismGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("AntiprismGraph",
        "AntiprismGraph[n] gives the n-antiprism: two offset n-cycles with cross "
        "edges (4-regular, 4n edges; n=3 is the octahedron).");

    symtab_add_builtin("PrismGraph", builtin_prism_graph);
    symtab_get_def("PrismGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PrismGraph",
        "PrismGraph[n] gives the n-gonal prism: two n-cycles joined by rungs "
        "(3-regular, 3n edges; n=4 is the cube).");

    symtab_add_builtin("SunletGraph", builtin_sunlet_graph);
    symtab_get_def("SunletGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("SunletGraph",
        "SunletGraph[n] gives the n-sunlet: a cycle C_n with a pendant vertex "
        "attached to each cycle vertex (2n vertices, 2n edges).");

    symtab_add_builtin("HelmGraph", builtin_helm_graph);
    symtab_get_def("HelmGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("HelmGraph",
        "HelmGraph[n] gives the helm graph: a wheel on an n-cycle rim with a "
        "pendant vertex attached to each rim vertex (2n+1 vertices, 3n edges).");

    symtab_add_builtin("GearGraph", builtin_gear_graph);
    symtab_get_def("GearGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GearGraph",
        "GearGraph[n] gives the gear graph: a hub joined to alternate vertices of "
        "a 2n-cycle (2n+1 vertices, 3n edges; bipartite).");

    symtab_add_builtin("EdgeBetweennessCentrality", builtin_edge_betweenness_centrality);
    symtab_get_def("EdgeBetweennessCentrality")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("EdgeBetweennessCentrality",
        "EdgeBetweennessCentrality[g] gives, for each edge, the number of "
        "shortest paths through it (exact rationals when paths tie).");

    symtab_add_builtin("DodecahedralGraph", builtin_dodecahedral_graph);
    symtab_get_def("DodecahedralGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DodecahedralGraph",
        "DodecahedralGraph[] gives the dodecahedron graph (GP(10,2)): 20 "
        "vertices, 30 edges, 3-regular.");

    symtab_add_builtin("IcosahedralGraph", builtin_icosahedral_graph);
    symtab_get_def("IcosahedralGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("IcosahedralGraph",
        "IcosahedralGraph[] gives the icosahedron graph: 12 vertices, 30 edges, "
        "5-regular.");

    /* ---- Phase 6: visualization ------------------------------------------ */
    symtab_add_builtin("GraphPlot", builtin_graph_plot);
    symtab_get_def("GraphPlot")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GraphPlot",
        "GraphPlot[g] gives a Graphics object drawing the graph g. Options: "
        "GraphLayout->\"name\" (e.g. \"SpringElectricalEmbedding\", "
        "\"CircularEmbedding\", \"SpiralEmbedding\", \"LinearEmbedding\", "
        "\"GridEmbedding\", \"RadialEmbedding\", \"LayeredEmbedding\", "
        "\"BipartiteEmbedding\", \"StarEmbedding\", \"RandomEmbedding\"), "
        "VertexStyle->color, EdgeStyle->color, VertexSize->r, "
        "VertexLabels->None. A bare Graph also auto-renders with these "
        "defaults.");

    symtab_add_builtin("Graph3D", builtin_graph3d);
    symtab_get_def("Graph3D")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Graph3D",
        "Graph3D[v, e] / Graph3D[e] builds a graph like Graph but displays it as "
        "a 3D node-link diagram (force-directed layout in a cube). Same edge "
        "sugar (u->v, u<->v) and simple-graph rules as Graph.");

    symtab_add_builtin("HighlightGraph", builtin_highlight_graph);
    symtab_get_def("HighlightGraph")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("HighlightGraph",
        "HighlightGraph[g, parts] draws g with the given vertices and/or edges "
        "emphasized (accent color, rest dimmed). Each part may be a vertex, an "
        "edge (u<->v / u->v), or a list of vertices treated as a path "
        "(highlighting its vertices and joining edges). Returns a Graphics.");
}

/* gmet_init.c - registration of the graph metrics builtins.
 *
 * Called once from the end of graph_init(). GraphDistance and CompleteGraph
 * are RE-registered here with wrappers that add the new argument forms
 * (GraphDistance[g, s], CompleteGraph[{n1, ...}]) and delegate every other
 * form to the original builtins, so neither shortestpath.c nor generators.c
 * changes.
 */

#include "graph_metrics.h"
#include "graph.h"
#include "symtab.h"
#include "attr.h"

static void reg(const char* name, BuiltinFunc f, const char* doc) {
    symtab_add_builtin(name, f);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(name, doc);
}

void graph_metrics_init(void) {
    /* ---- distances ------------------------------------------------------- */
    reg("GraphDistanceMatrix", builtin_graph_distance_matrix,
        "GraphDistanceMatrix[g] gives the matrix of shortest-path distances "
        "between all pairs of vertices of g, rows and columns in VertexList "
        "order; entry {i, j} is the distance from vertex i to vertex j "
        "(Infinity if unreachable). GraphDistanceMatrix[g, d] keeps only "
        "distances at most d. Integers for unweighted graphs, machine reals "
        "when g has EdgeWeight (edge weights as lengths).");
    reg("GraphDistance", builtin_gmet_graph_distance,
        "GraphDistance[g, s, t] gives the length of a shortest path from s to t "
        "(Infinity if unreachable). GraphDistance[g, s] gives the list of "
        "distances from s to every vertex, in VertexList order. Edge weights "
        "are used as lengths when g has EdgeWeight.");
    reg("VertexEccentricity", builtin_vertex_eccentricity,
        "VertexEccentricity[g, v] gives the largest distance from v to any "
        "vertex reachable from v. With EdgeWeight the weights are lengths and, "
        "as in Wolfram, a vertex v cannot reach makes it Infinity.");
    reg("GraphDiameter", builtin_graph_diameter,
        "GraphDiameter[g] gives the greatest distance between two vertices of "
        "g: the maximum vertex eccentricity. Infinity unless g is (strongly) "
        "connected.");
    reg("GraphRadius", builtin_graph_radius,
        "GraphRadius[g] gives the minimum vertex eccentricity of g. Infinity "
        "unless g is (strongly) connected.");
    reg("GraphCenter", builtin_graph_center,
        "GraphCenter[g] gives the vertices of g with minimum eccentricity; {} "
        "unless g is (strongly) connected.");
    reg("GraphPeriphery", builtin_graph_periphery,
        "GraphPeriphery[g] gives the vertices of g with maximum eccentricity; "
        "{} unless g is (strongly) connected.");
    reg("MeanGraphDistance", builtin_mean_graph_distance,
        "MeanGraphDistance[g] gives the mean distance over all ordered pairs of "
        "distinct vertices of g: exact for unweighted graphs, a machine real "
        "for weighted ones, Infinity unless g is (strongly) connected.");
    reg("GraphDensity", builtin_graph_density,
        "GraphDensity[g] gives the number of edges of g divided by the number "
        "of possible edges: (directed edges + 2 undirected edges)/(n(n-1)).");
    reg("KirchhoffMatrix", builtin_kirchhoff_matrix,
        "KirchhoffMatrix[g] gives the Kirchhoff (Laplacian) matrix D - A of g, "
        "with D the diagonal matrix of vertex degrees (incident edges) and A "
        "the adjacency matrix. Dense (Wolfram returns a SparseArray); weights "
        "are ignored.");

    /* ---- centralities ---------------------------------------------------- */
    reg("DegreeCentrality", builtin_degree_centrality,
        "DegreeCentrality[g] gives the list of vertex degrees of g; "
        "DegreeCentrality[g, \"In\"] and [g, \"Out\"] count incoming and "
        "outgoing edges. In a mixed graph an undirected edge counts as one edge "
        "in each direction.");
    reg("ClosenessCentrality", builtin_closeness_centrality,
        "ClosenessCentrality[g] gives, for each vertex v, r/s where r is the "
        "number of vertices reachable from v and s the sum of their distances "
        "from v (0 if v reaches none). Uses EdgeWeight as lengths.");
    reg("EccentricityCentrality", builtin_eccentricity_centrality,
        "EccentricityCentrality[g] gives 1/VertexEccentricity[g, v] for each "
        "vertex v (0 when the eccentricity is 0). Uses EdgeWeight as lengths.");
    reg("BetweennessCentrality", builtin_betweenness_centrality,
        "BetweennessCentrality[g] gives for each vertex v the sum over pairs of "
        "other vertices s, t of the fraction of shortest s-t paths passing "
        "through v (unordered pairs for undirected graphs, ordered pairs for "
        "directed ones). Brandes' algorithm; edge weights are ignored, as in "
        "Wolfram. Mixed graphs are not supported.");
    reg("EdgeBetweennessCentrality", builtin_edge_betweenness_centrality,
        "EdgeBetweennessCentrality[g] gives for each edge, in EdgeList order, "
        "the sum over ordered pairs of vertices s, t of the fraction of "
        "shortest s-t paths using that edge. Uses EdgeWeight as lengths.");
    reg("PageRankCentrality", builtin_pagerank_centrality,
        "PageRankCentrality[g, a] gives the PageRank of each vertex with "
        "damping factor a (default 0.85): x = a P^T x + (1-a)/n, where a "
        "vertex with no outgoing edge links to every vertex; the entries sum "
        "to 1. Edge weights are ignored.");
    reg("EigenvectorCentrality", builtin_eigenvector_centrality,
        "EigenvectorCentrality[g] gives the eigenvector centrality of each "
        "vertex, computed per strongly connected component (each component of "
        "k > 1 vertices weighted by k - 1, total 1; isolated vertices 0). "
        "EigenvectorCentrality[g, \"In\"] (default) uses incoming edges, "
        "[g, \"Out\"] outgoing ones. Edge weights are ignored.");
    reg("KatzCentrality", builtin_katz_centrality,
        "KatzCentrality[g, a] gives the Katz centrality x = a A^T x + 1 of "
        "each vertex; KatzCentrality[g, a, b] uses b (a number or a list) in "
        "place of 1. Edge weights are ignored.");
    reg("HITSCentrality", builtin_hits_centrality,
        "HITSCentrality[g] gives {h, a}: h is the principal eigenvector of "
        "A^T A (per block of vertices sharing in-neighbours, each block of "
        "k > 1 vertices weighted k - 1, total 1) and a = A h, with A the "
        "adjacency matrix (Wolfram's convention). Edge weights are ignored.");

    /* ---- clustering ------------------------------------------------------ */
    reg("GraphTriangleCount", builtin_graph_triangle_count,
        "GraphTriangleCount[g] gives the number of triangles of an undirected "
        "graph g, or of directed 3-cycles of a directed graph. O(m^1.5) "
        "degree-ordered triangle listing.");
    reg("LocalClusteringCoefficient", builtin_local_clustering_coefficient,
        "LocalClusteringCoefficient[g] gives for each vertex the fraction of "
        "pairs of its neighbours that are adjacent (exact); "
        "LocalClusteringCoefficient[g, v] gives it for v. For directed graphs, "
        "directed 3-cycles through v over in/out neighbour pairs.");
    reg("GlobalClusteringCoefficient", builtin_global_clustering_coefficient,
        "GlobalClusteringCoefficient[g] gives 3 x (number of triangles) / "
        "(number of connected triples) of g, exactly.");
    reg("MeanClusteringCoefficient", builtin_mean_clustering_coefficient,
        "MeanClusteringCoefficient[g] gives the mean of the local clustering "
        "coefficients of g, exactly.");

    /* ---- generators ------------------------------------------------------ */
    reg("WheelGraph", builtin_wheel_graph,
        "WheelGraph[n] gives the wheel graph with n vertices: vertex 1 joined "
        "to every vertex of the cycle 2, ..., n.");
    reg("HypercubeGraph", builtin_hypercube_graph,
        "HypercubeGraph[n] gives the n-dimensional hypercube graph on 2^n "
        "vertices.");
    reg("GridGraph", builtin_grid_graph,
        "GridGraph[{m, n}] gives the m x n grid graph; GridGraph[{n1, ..., nk}] "
        "the k-dimensional grid. The first coordinate varies fastest in the "
        "vertex numbering.");
    reg("KaryTree", builtin_kary_tree,
        "KaryTree[n] gives a binary tree with n vertices; KaryTree[n, k] a "
        "k-ary tree with n vertices, vertices numbered in breadth-first order.");
    reg("CompleteKaryTree", builtin_complete_kary_tree,
        "CompleteKaryTree[n] gives the complete binary tree with n levels; "
        "CompleteKaryTree[n, k] the complete k-ary tree with n levels.");
    reg("CirculantGraph", builtin_circulant_graph,
        "CirculantGraph[n, j] gives the circulant graph on n vertices with "
        "i joined to i+j and i-j (mod n); CirculantGraph[n, {j1, j2, ...}] "
        "uses every jump ji.");
    reg("PetersenGraph", builtin_petersen_graph,
        "PetersenGraph[] gives the Petersen graph; PetersenGraph[n, k] the "
        "generalized Petersen graph: inner vertices 1..n with i joined to i+k "
        "(mod n), outer cycle n+1..2n, and spokes i to n+i.");
    reg("TuranGraph", builtin_turan_graph,
        "TuranGraph[n, k] gives the Turan graph: the complete k-partite graph "
        "on n vertices with parts as equal as possible.");
    reg("HararyGraph", builtin_harary_graph,
        "HararyGraph[k, n] gives the Harary graph: a k-connected graph on n "
        "vertices with the minimum number of edges (k >= 2, n > k).");
    reg("CompleteGraph", builtin_gmet_complete_graph,
        "CompleteGraph[n] gives the complete graph on n vertices. "
        "CompleteGraph[{n1, n2, ...}] gives the complete multipartite graph "
        "with parts of sizes n1, n2, ....");
}

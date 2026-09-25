# WeightedAdjacencyMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WeightedAdjacencyMatrix[g] gives the adjacency matrix of g with each entry the corresponding edge's weight (0 where there is no edge). Equal to AdjacencyMatrix[g] when g has no EdgeWeight.`**

## Examples

_No verified examples yet for this function._

## Algorithm

wtadjmat.c - WeightedAdjacencyMatrix[g]: dense adjacency matrix filled with per-edge weights instead of a literal 1.

Same algorithm as AdjacencyMatrix (adjmat.c): a DirectedEdge[a,b] sets M[a][b] = weight(a,b); an UndirectedEdge sets both M[a][b] and M[b][a]. Any entry with no edge is 0. For a graph with no EdgeWeight, every weight defaults to 1 (graph_resolve_edge_weights), so WeightedAdjacencyMatrix[g] == AdjacencyMatrix[g] exactly for an unweighted g.

Memory (SPEC section 4): returns a freshly-allocated matrix; frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

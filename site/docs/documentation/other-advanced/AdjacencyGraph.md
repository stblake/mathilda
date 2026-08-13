# AdjacencyGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AdjacencyGraph[m] builds a graph on vertices 1..n from a 0/1 adjacency matrix m (undirected if m is symmetric, else directed).`**

## Examples

_No verified examples yet for this function._

## Algorithm

adjgraph.c - AdjacencyGraph[m]: build a graph from a 0/1 adjacency matrix.

The inverse of AdjacencyMatrix. Vertices are the integers 1..n. A symmetric matrix yields an undirected graph (one UndirectedEdge per i<j with m[i][j]=1); an asymmetric matrix yields a directed graph (a DirectedEdge for each off- diagonal m[i][j]=1). Diagonal entries (self-loops) are ignored. The result is returned as a Graph[...] expression and canonicalized/validated by the evaluator (builtin_graph).

Round-trips with AdjacencyMatrix when the source graph's vertices are 1..n.

Memory (SPEC section 4): returns a freshly-allocated Graph; frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

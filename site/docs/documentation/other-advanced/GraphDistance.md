# GraphDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphDistance[g,s,t] gives the length of a shortest path from s to t (Infinity if unreachable).`**

## Examples

_No verified examples yet for this function._

## Algorithm

shortestpath.c - FindShortestPath[g,s,t] and GraphDistance[g,s,t].

Unweighted breadth-first search over the successor adjacency (out[]): for a directed graph this follows edge direction; for an undirected graph out[] is symmetric, so it is an ordinary shortest path. Wolfram's naming split is kept: FindShortestPath returns the vertex path, GraphDistance the length.

Unreachable target: FindShortestPath -> {} (empty list), GraphDistance -> Infinity.

Memory (SPEC section 4): returns freshly-allocated results; frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

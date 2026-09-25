# FindShortestPath

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindShortestPath[g,s,t] gives a shortest path from s to t as a list of vertices ({} if none).`**

## Examples

_No verified examples yet for this function._

## Algorithm

shortestpath.c - FindShortestPath[g,s,t] and GraphDistance[g,s,t].

Two algorithms, dispatched on graph_weights_usable(g):

```text
  - Unweighted (default): breadth-first search over the successor adjacency
    (GraphAdj.out[]): for a directed graph this follows edge direction; for
    an undirected graph out[] is symmetric, so it is an ordinary shortest
    path.
  - Weighted (g carries a non-negative-numeric EdgeWeight): Dijkstra over a
    local, call-scoped weighted adjacency (WAdj below) -- NOT over
    GraphAdj, which has no weight storage and is shared by 5 other
    builtins (ConnectedComponents, WeaklyConnectedComponents,
    FindSpanningTree, ConnectedGraphQ, VertexConnectivity); widening it
    would risk the exact class of shared-choke-point defect a
    plan-reviewer pass caught during the EdgeWeight ticket. Falls back to
    BFS for a symbolic or negative weight rather than erroring.
```

Wolfram's naming split is kept: FindShortestPath returns the vertex path, GraphDistance the length/total weight.

Unreachable target: FindShortestPath -> {} (empty list), GraphDistance -> Infinity.

Memory (SPEC section 4): returns freshly-allocated results; frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

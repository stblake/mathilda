# GraphDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphDistance[g, s, t] gives the length of a shortest path from s to t (Infinity if unreachable). GraphDistance[g, s] gives the list of distances from s to every vertex, in VertexList order. Edge weights are used as lengths when g has EdgeWeight.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= GraphDistance[Graph[{1->2, 2->3, 3->1, 3->4}], 1]
Out[1]= {0, 1, 2, 3}

In[2]:= GraphDistance[Graph[{1->2, 2->3, 3->1, 3->4}], 4, 1]
Out[2]= Infinity
```

### Options (3)

```mathematica
In[3]:= GraphDistance[Graph[{1,2,3},{1<->2,2<->3,1<->3}, EdgeWeight->{1,1,5}], 1]
Out[3]= {0, 1.0, 2.0}

In[4]:= GraphDistance[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{5,7}], 1, 3]
Out[4]= 12.0

In[5]:= GraphDistance[Graph[{1,2},{1<->2}, EdgeWeight->{a}], 1]
Out[5]= GraphDistance[Graph[<2 vertices, 1 edge>], 1]
```

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

- *(w)* weight-aware. `GraphDistance[g, s]` is the single-source form added
  by the graph-metrics module (`src/graph/gmet_*.c`); `GraphDistance[g, s, t]`
  is unchanged. `GraphDistance` is re-registered by a wrapper that delegates
  its pre-existing form (`[g, s, t]`) to the original builtin.
- `[g, s, t]` uses the same weight-aware dispatch as `FindShortestPath`
  (Dijkstra when every weight is non-negative and numeric, BFS otherwise) and
  on a weighted graph returns a machine real, as the Wolfram Language does
  (weights `{5, 7}` give `12.`), identical to `GraphDistance[g, s]` and
  `GraphDistanceMatrix`.
- Single-source form, weighted: machine reals, except the source's own entry,
  which is an exact `0` as in Wolfram. A symbolic, complex or negative weight
  leaves the single-source form unevaluated.
- Unweighted single-source results are Integers, packed when every vertex is
  reachable. Edge direction is followed; an undirected edge is usable both
  ways.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/), [FindShortestPath](../../graphs/FindShortestPath/), [GraphDistanceMatrix](../../graphs/GraphDistanceMatrix/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

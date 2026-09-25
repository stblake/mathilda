# FindShortestPath

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindShortestPath[g,s,t] gives a shortest path from s to t as a list of vertices ({} if none).`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= FindShortestPath[Graph[{1,2,3,4},{1->2,2->3,3->4}], 1, 4]
Out[1]= {1, 2, 3, 4}

In[2]:= FindShortestPath[Graph[{1,2,3,4},{1->2,2->3,3->4}], 4, 1]
Out[2]= {}

In[3]:= FindShortestPath[CycleGraph[6], 1, 4]
Out[3]= {1, 2, 3, 4}
```

### Scope (1)

```mathematica
In[4]:= GraphDistance[Graph[{1,2,3,4},{1->2,2->3,3->4}], 4, 1]
Out[4]= Infinity
```

### Options (3)

```mathematica
In[5]:= FindShortestPath[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,10}], 1, 4]
Out[5]= {1, 2, 3, 4}

In[6]:= FindShortestPath[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,-10}], 1, 4]
Out[6]= {1, 4}

In[7]:= GraphDistance[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,10}], 1, 4]
Out[7]= 3.0
```

## Options & behaviour

### Weighted

the direct `1 -> 4` edge (weight 10) loses to the longer, cheaper
`1 -> 2 -> 3 -> 4` route (weight 3). A negative weight falls back to hop count.

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

- `Protected`. Shared by the search & computation heads (`FindShortestPath`,
  `GraphDistance`, `ConnectedComponents`, `WeaklyConnectedComponents`,
  `StronglyConnectedComponents`, `FindSpanningTree`, `ConnectedGraphQ`,
  `VertexConnectivity`): all build an integer-indexed adjacency on demand, and
  all but `FindShortestPath`/`GraphDistance` are unweighted.
- **Weight-aware**: if `g` carries an `EdgeWeight` and every weight is
  non-negative and numeric, uses Dijkstra (minimum total weight); otherwise
  (unweighted, a symbolic weight, or a negative weight present) uses unweighted
  BFS (minimum hop count), following edge direction for directed graphs either
  way.
- Weighted search is a plain `O(V^2)` Dijkstra (no priority queue — consistent
  with `VertexConnectivity`'s original small-graph exact-algorithm precedent),
  and it falls back to unweighted BFS rather than erroring whenever a weight
  isn't usable for it (not present, symbolic, or negative). There is no
  Bellman-Ford / negative-weight support.
- Unevaluated when `s` or `t` is not a vertex.
- Companion `GraphDistance[g, s, t]` gives the length/total weight of that path;
  `Infinity` if unreachable. Same weight-aware dispatch as `FindShortestPath`,
  and, on a weighted graph, returns a machine real as the Wolfram Language does
  (weights `{5, 7}` give `12.`), identical to `GraphDistance[g, s]` and
  `GraphDistanceMatrix`. (Until v0.186 it returned an exact
  `Integer`/`Rational`, which disagreed with both Mathematica and the
  single-source form.)

**Attributes:** `Protected`.

## References

**See also:** [GraphDistance](../../graphs/GraphDistance/), [ConnectedComponents](../../graphs/ConnectedComponents/), [WeaklyConnectedComponents](../../graphs/WeaklyConnectedComponents/), [StronglyConnectedComponents](../../graphs/StronglyConnectedComponents/), [FindSpanningTree](../../graphs/FindSpanningTree/), [ConnectedGraphQ](../../graphs/ConnectedGraphQ/), [VertexConnectivity](../../graphs/VertexConnectivity/), [EdgeWeight](../../graphs/EdgeWeight/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

# FindShortestPath

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindShortestPath[g,s,t] gives a shortest path from s to t as a list of vertices ({} if none).`**

## Examples (9)

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

### Applications (2)

Follows edge directions

```mathematica
In[8]:= FindShortestPath[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}], 1, 4]
Out[8]= {1, 2, 3, 4}
```

No path exists

```mathematica
In[9]:= FindShortestPath[Graph[{1 -> 2, 3 -> 4}], 1, 4]
Out[9]= {}
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

**Algorithm.** `builtin_find_shortest_path` returns a shortest `s`-`t` path as a list of
vertices. It dispatches on `graph_weights_usable(g)`: a graph carrying an `EdgeWeight` list in
which every weight is a non-negative, non-complex number runs **Dijkstra**; otherwise — an
unweighted graph, or one with any symbolic, negative or complex weight — it falls back to
**BFS** rather than erroring. Both follow edge direction on a directed graph and treat an
undirected edge as usable both ways, and both reconstruct the path by walking a `parent[]` array
back from `t`.

**Data structures.** BFS uses the shared CSR `GraphAdj` with an integer FIFO queue, a `dist[]`
array (`-1` = unreached) and `parent[]`. Dijkstra builds a *separate* call-scoped weighted
adjacency `WAdj` (via a `GraphVIdx` hash) rather than widening the shared `GraphAdj`, with a
`double` `dist[]` (`DBL_MAX` = unreached), a `done[]` flag array and `parent[]`; it selects the
next vertex by a linear array scan, not a heap.

**Complexity / limits.** BFS is `O(V + E)`; the array-scan Dijkstra is `O(V^2)`. When there is
no `s`-`t` path the result is `{}` (the empty list). A non-graph argument, or an `s`/`t` that is
not a vertex, returns unevaluated. (The companion `GraphDistance` returns the path *length*, or
`Infinity` when unreachable.)

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

- E. W. Dijkstra, *A note on two problems in connexion with graphs*, Numer. Math. **1** (1959) 269-271.
- Source: [`src/graph/shortestpath.c`](https://github.com/stblake/mathilda/blob/main/src/graph/shortestpath.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

The result is the path from `s` to `t` as a list of vertices; it is `{}` when no path exists.
On a directed graph the path respects edge directions.

A graph that carries non-negative numeric `EdgeWeight`s is traversed with Dijkstra's algorithm,
so the path minimises total weight; an unweighted graph (or one with symbolic or negative
weights) falls back to a breadth-first search that minimises the hop count. `GraphDistance`
gives the corresponding path length.

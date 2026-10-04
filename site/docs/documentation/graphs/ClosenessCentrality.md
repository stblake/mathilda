# ClosenessCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ClosenessCentrality[g] gives, for each vertex v, r/s where r is the number of vertices reachable from v and s the sum of their distances from v (0 if v reaches none). Uses EdgeWeight as lengths.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ClosenessCentrality[PathGraph[4]]
Out[1]= {0.5, 0.75, 0.75, 0.5}

In[2]:= ClosenessCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= {0.5, 0.6, 0.75, 0.0}

In[3]:= ClosenessCentrality[Graph[{1,2,3},{}]]
Out[3]= {0.0, 0.0, 0.0}
```

### Options (1)

```mathematica
In[4]:= ClosenessCentrality[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[4]= {0.5, 0.5, 0.0}
```

### Applications (2)

The middle vertex is closest to the rest

```mathematica
In[5]:= ClosenessCentrality[Graph[{1 <-> 2, 2 <-> 3}]]
Out[5]= {0.666667, 1.0, 0.666667}
```

Directed: only forward-reachable vertices count

```mathematica
In[6]:= ClosenessCentrality[Graph[{1 -> 2, 2 -> 3, 1 -> 3, 3 -> 4}]]
Out[6]= {0.75, 0.666667, 1.0, 0.0}
```

## Implementation notes

**Algorithm.** `builtin_closeness_centrality` reduces from a cached per-source distance
summary, `gmet_dist_summary(g)`, which runs a multi-source traversal — BFS when the graph is
unweighted, Dijkstra when it carries usable `EdgeWeight` lengths — from every vertex. For each
vertex `v` it records `reach[v]`, the number of vertices `v` can reach, and `sum[v]`, the total
distance to those vertices. The centrality is then `v = reach[v] / sum[v]` (and `0` when `v`
reaches nothing). Dividing by the sum over only the *reachable* set keeps the value finite on a
disconnected graph, which is the practical difference from the textbook `(n-1)/sum` definition.

**Data structures.** The shared summary holds the `reach[]` and `sum[]` arrays (and the
eccentricity used by `GraphCenter`/`GraphRadius`), memoised on the graph node so repeated metric
queries on the same graph are `O(1)`. Distances follow edge direction on a directed graph and an
undirected edge is usable both ways. The output is a packed machine-real vector.

**Complexity / limits.** The summary costs `O(V(V+E))` unweighted and `O(V(E + V log V))`
weighted; the reduction to the centrality vector is `O(V)`. Single argument only; non-graph
input, or a summary that fails to build, returns unevaluated.

- *(w)* weight-aware; machine reals, packed.
- For a vertex `v`: `r/s`, with `r` the number of vertices reachable from `v`
  and `s` the sum of their distances; 0 if none are reachable.
- Reduced from the cached MS-BFS per-source summary (weighted: Dijkstra per
  source); see `GraphDistanceMatrix`.

**Attributes:** `Protected`.

## References

**See also:** [GraphDistanceMatrix](../../graphs/GraphDistanceMatrix/)

- Source: [`src/graph/gmet_centrality.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_centrality.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

For each vertex `v` the value is `r/s`, where `r` is the number of vertices reachable from `v`
and `s` is the sum of their distances from `v`; a vertex that reaches nothing gets `0`.
Dividing by the reachable set (rather than by `n-1`) keeps the measure finite on a disconnected
or directed graph.

`EdgeWeight` is used as edge length when present, so on a weighted graph the distances — and
hence the centralities — reflect the weights rather than hop counts.

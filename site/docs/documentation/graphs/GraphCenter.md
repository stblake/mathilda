# GraphCenter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphCenter[g] gives the vertices of g with minimum eccentricity; {} unless g is (strongly) connected.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= GraphCenter[PathGraph[5]]
Out[1]= {3}

In[2]:= GraphPeriphery[PathGraph[5]]
Out[2]= {1, 5}

In[3]:= GraphCenter[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= {}

In[4]:= GraphPeriphery[Graph[{1<->2,3<->4}]]
Out[4]= {}
```

### Options (1)

```mathematica
In[5]:= GraphCenter[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[5]= {1}
```

### Applications (1)

The one middle vertex of an odd path

```mathematica
In[6]:= GraphCenter[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5}]]
Out[6]= {3}
```

## Implementation notes

**Algorithm.** `builtin_graph_center` gives the vertices whose eccentricity equals the graph
radius (the minimum eccentricity). Through the shared `extremal_compute`, it reads the cached
per-source distance summary (`gmet_dist_summary`: multi-source BFS unweighted, Dijkstra
weighted), finds the minimum eccentricity `target`, and returns every vertex whose eccentricity
matches it (within a relative tolerance `ECC_TIE_TOL = 1e-12` on the weighted path). An
`O(V+E)` strong-connectivity test short-circuits first: for an unweighted graph that is not
(strongly) connected the centre is empty.

**Data structures.** The memoised distance summary's `ecc[]` array, shared with `GraphRadius`
and `MeanGraphDistance`; `gmet_vertex_subset` lifts the flagged indices back to the vertex
expressions in `VertexList` order.

**Complexity / limits.** `O(V(V+E))` unweighted, `O(V(E + V log V))` weighted. The empty graph
and a non-(strongly-)connected unweighted graph both return `{}`; a non-graph argument returns
unevaluated.

- *(w)* weight-aware.
- Unweighted and not strongly connected (not connected, if undirected): `{}`.
- Weighted: taken over the weighted eccentricities, so a weighted digraph that
  is not strongly connected can have a non-empty center
  (*reverse-engineered*).
- Reduced from the cached MS-BFS per-source summary (see `GraphDistanceMatrix`).

**Attributes:** `Protected`.

## References

**See also:** [GraphPeriphery](../../graphs/GraphPeriphery/), [GraphDistanceMatrix](../../graphs/GraphDistanceMatrix/)

- Source: [`src/graph/gmet_distance.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_distance.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The centre is the set of vertices of minimum eccentricity — those whose greatest distance to any
other vertex equals the graph radius. It is returned as a list of vertices in `VertexList`
order.

It is `{}` unless the graph is connected (strongly connected, for a directed graph). The
matching minimum eccentricity value is `GraphRadius`.

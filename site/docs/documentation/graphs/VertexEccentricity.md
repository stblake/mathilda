# VertexEccentricity

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexEccentricity[g, v] gives the largest distance from v to any vertex reachable from v. With EdgeWeight the weights are lengths and, as in Wolfram, a vertex v cannot reach makes it Infinity.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= VertexEccentricity[PathGraph[5], 3]
Out[1]= 2

In[2]:= VertexEccentricity[Graph[{1->2, 2->3, 3->1, 3->4}], 1]
Out[2]= 3

In[3]:= VertexEccentricity[Graph[{1<->2,3<->4}], 1]
Out[3]= 1
```

### Options (2)

```mathematica
In[4]:= VertexEccentricity[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}], 1]
Out[4]= 3.0

In[5]:= VertexEccentricity[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}], 2]
Out[5]= Infinity
```

### Applications (3)

An end is four hops from the far end

```mathematica
In[6]:= VertexEccentricity[PathGraph[{1, 2, 3, 4, 5}], 1]
Out[6]= 4
```

The middle vertex is closer to everything

```mathematica
In[7]:= VertexEccentricity[PathGraph[{1, 2, 3, 4, 5}], 3]
Out[7]= 2
```

The farthest vertex on a 5-cycle is two hops away

```mathematica
In[8]:= VertexEccentricity[CycleGraph[5], 1]
Out[8]= 2
```

## Implementation notes

**Algorithm.** `builtin_vertex_eccentricity` gives `VertexEccentricity[g, v]`, the
largest distance from `v` to any other vertex. It runs a single-source distance
computation from `v` (`single_source` — BFS for an unweighted graph, Dijkstra
when `g` carries usable `EdgeWeight`) and takes the maximum. For an **unweighted**
graph the maximum is over the vertices `v` can reach; for a **weighted** graph
Wolfram measures over *all* vertices, so a vertex `v` cannot reach makes the
eccentricity `Infinity`. The result is an `Integer` for an unweighted graph and a
machine `Real` otherwise.

**Data structures.** The CSR distance machinery of `src/graph/gmet_distance.c`:
an `int64 di[]` array of hop distances and a `double dw[]` of weighted distances,
filled by the single-source routine from the memo's integer endpoints and the
resolved edge weights (`gmet_edge_weights`).

**Complexity / limits.** `O(V + E)` unweighted, `O(E + V log V)` weighted, per
call. The graph-wide maxima/minima of eccentricity are `GraphDiameter`,
`GraphRadius`, `GraphCenter`, and `GraphPeriphery`, which compute all
eccentricities at once and cache them. A `v` that is not a vertex leaves the call
unevaluated.

- *(w)* weight-aware (machine reals when weighted; symbolic, complex or
  negative weights leave it unevaluated).
- Unweighted: the largest distance to a vertex that `v` reaches, so it is
  finite on disconnected graphs.
- Weighted: `Infinity` if `v` does not reach every vertex (Wolfram's weighted
  rule).
- Computed from the cached per-source MS-BFS summary (see
  `GraphDistanceMatrix`).

**Attributes:** `Protected`.

## References

**See also:** [GraphDistanceMatrix](../../graphs/GraphDistanceMatrix/)

- Source: [`src/graph/gmet_distance.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_distance.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The eccentricity of a vertex is the greatest shortest-path distance from it to
any other vertex. The minimum eccentricity over the graph is the radius, the
maximum is the diameter, and the vertices achieving them form the center and
periphery.

For an unweighted graph the answer is an integer. With `EdgeWeight` the weights
are lengths, and a vertex that cannot reach every other has eccentricity
`Infinity`.

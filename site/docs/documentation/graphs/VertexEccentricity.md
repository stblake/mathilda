# VertexEccentricity

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexEccentricity[g, v] gives the largest distance from v to any vertex reachable from v. With EdgeWeight the weights are lengths and, as in Wolfram, a vertex v cannot reach makes it Infinity.`**

## Examples (5)

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

## Implementation notes

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

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

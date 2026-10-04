# LocalClusteringCoefficient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LocalClusteringCoefficient[g] gives for each vertex the fraction of pairs of its neighbours that are adjacent (exact); LocalClusteringCoefficient[g, v] gives it for v. For directed graphs, directed 3-cycles through v over in/out neighbour pairs.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= LocalClusteringCoefficient[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[1]= {1, 1, 1/3, 0}

In[2]:= LocalClusteringCoefficient[Graph[{1<->2,2<->3,3<->1,3<->4}], 3]
Out[2]= 1/3

In[3]:= LocalClusteringCoefficient[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= {1, 1, 1/2, 0}

In[4]:= LocalClusteringCoefficient[Graph[{1<->2, 2->3}]]
Out[4]= LocalClusteringCoefficient[Graph[<3 vertices, 2 edges>]]
```

### Options (1)

```mathematica
In[5]:= LocalClusteringCoefficient[Graph[{1,2,3,4},{1<->2,2<->3,3<->1,3<->4}, EdgeWeight->{1,2,3,4}]]
Out[5]= {1, 1, 1/3, 0}
```

### Applications (3)

Every neighbourhood is complete, so all 1

```mathematica
In[6]:= LocalClusteringCoefficient[CompleteGraph[4]]
Out[6]= {1, 1, 1, 1}
```

No neighbour pair is adjacent

```mathematica
In[7]:= LocalClusteringCoefficient[CycleGraph[5]]
Out[7]= {0, 0, 0, 0, 0}
```

The hub: 4 of 6 neighbour pairs are joined

```mathematica
In[8]:= LocalClusteringCoefficient[WheelGraph[5], 1]
Out[8]= 2/3
```

## Implementation notes

**Algorithm.** `builtin_local_clustering_coefficient` gives, per vertex, the
fraction of pairs of its neighbours that are adjacent: `t(v) / C(d(v), 2)`
(`0` when `d(v) < 2`), where `t(v)` is the number of triangles through `v`. For a
directed graph it counts directed 3-cycles through `v` over the `in(v) out(v) -
r(v)` in/out neighbour pairs. `LocalClusteringCoefficient[g]` returns the whole
list in `VertexList` order; `[g, v]` the single value. The per-vertex triangle
counts come from the same degree-ordered triangle listing as
`GraphTriangleCount`, run with `per_vertex = 1`.

**Data structures.** The CSR oriented adjacency and thread-team triangle listing
of `tri_compute` (`TriInfo.t[v]`, `TriInfo.den[v]`). Because degree sequences
repeat, most coefficients repeat, so the whole-graph form builds each distinct
`(t, den)` value once and shares it by reference through a small open-addressed
hash; the result is cached on the graph node.

**Complexity / limits.** `O(m^1.5)` for the listing, `O(V)` to assemble. Values
are exact `Rational`s (or `0`). `EdgeWeight` is ignored; a mixed graph is left
unevaluated.

- All clustering heads are exact, ignore weights, and leave mixed graphs
  unevaluated (as in Wolfram).
- Undirected: `t(v)/C(d(v), 2)`, with `t(v)` the triangles at `v` and `d(v)`
  its degree.
- Directed (*reverse-engineered*): triangles are directed 3-cycles; local
  coefficient `c(v)/(in(v) out(v) - r(v))`, with `c(v)` the directed 3-cycles
  through `v` and `r(v)` the reciprocally linked neighbours.
- Triangles are found by the oriented triangle-listing algorithm of
  `GraphTriangleCount`.

**Attributes:** `Protected`.

## References

**See also:** [GraphTriangleCount](../../graphs/GraphTriangleCount/)

- N. Chiba and T. Nishizeki, *Arboricity and subgraph listing algorithms*, SIAM J. Comput. **14** (1985) 210-223.
- Source: [`src/graph/gmet_cluster.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_cluster.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The local clustering coefficient of a vertex is how close its neighbourhood is
to a clique: the number of edges among its neighbours divided by the number of
possible such edges, `C(d, 2)`. A vertex of degree below 2 gets `0`.

`LocalClusteringCoefficient[g]` returns one exact `Rational` per vertex in
`VertexList` order; `LocalClusteringCoefficient[g, v]` gives the single value.
The mean of these numbers is `MeanClusteringCoefficient`.

# LocalClusteringCoefficient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LocalClusteringCoefficient[g] gives for each vertex the fraction of pairs of its neighbours that are adjacent (exact); LocalClusteringCoefficient[g, v] gives it for v. For directed graphs, directed 3-cycles through v over in/out neighbour pairs.`**

## Examples (5)

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

## Implementation notes

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

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

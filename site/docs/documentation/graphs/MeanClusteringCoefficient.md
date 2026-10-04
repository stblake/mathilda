# MeanClusteringCoefficient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MeanClusteringCoefficient[g] gives the mean of the local clustering coefficients of g, exactly.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= MeanClusteringCoefficient[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[1]= 7/12

In[2]:= MeanClusteringCoefficient[CompleteGraph[4]]
Out[2]= 1

In[3]:= MeanClusteringCoefficient[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= 5/8

In[4]:= MeanClusteringCoefficient[Graph[{1<->2, 2->3}]]
Out[4]= MeanClusteringCoefficient[Graph[<3 vertices, 2 edges>]]
```

### Applications (5)

Every neighbourhood is fully connected

```mathematica
In[5]:= MeanClusteringCoefficient[CompleteGraph[4]]
Out[5]= 1
```

A cycle of five has no triangles

```mathematica
In[6]:= MeanClusteringCoefficient[CycleGraph[5]]
Out[6]= 0
```

A triangle with a pendant vertex

```mathematica
In[7]:= MeanClusteringCoefficient[Graph[{1, 2, 3, 4}, {UndirectedEdge[1, 2], UndirectedEdge[2, 3], UndirectedEdge[1, 3], UndirectedEdge[3, 4]}]]
Out[7]= 7/12
```

The leaves have degree 1, so their local coefficient is 0

```mathematica
In[8]:= MeanClusteringCoefficient[StarGraph[5]]
Out[8]= 0
```

A directed 3-cycle counts as a triangle

```mathematica
In[9]:= MeanClusteringCoefficient[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3], DirectedEdge[3, 1]}]]
Out[9]= 1
```

## Implementation notes

**Algorithm.** `builtin_mean_clustering_coefficient` averages the local clustering coefficients. For an undirected graph the local value is `t(v) / C(d(v), 2)`, with `t(v)` the triangles through `v` and `d(v)` its degree, and `0` when `d(v) < 2`. Triangles are listed over the underlying simple graph with the degree-ordered orientation (each edge points from lower to higher `(degree, index)` rank). For a directed graph a triangle is a directed 3-cycle, and the local denominator is `in(v) out(v) - r(v)`, where `r(v)` counts neighbours joined in both directions. The directed convention was reverse-engineered from Mathematica 15. The mean is exact: numerators are grouped by denominator and summed as GMP rationals (`mpq_t`), then divided by `n`. The empty graph gives `0`.

**Data structures.** `tri_compute` fills per-vertex `int64_t` arrays `t` and `den`. Each oriented edge carries a two-bit record of which arcs exist, so the directed 3-cycle test is two bit-ANDs per triangle. Rows are processed on the thread team with per-thread mark arrays. The result is an Integer or Rational `Expr`.

**Complexity / limits.** `O(m^1.5)` total, since every oriented out-degree is bounded by `O(sqrt m)`. Weights are ignored, and a mixed graph (directed and undirected edges together) is left unevaluated.

- Exact, weights ignored, mixed graphs unevaluated (as in Wolfram).
- The arithmetic mean of `LocalClusteringCoefficient[g]` (undirected or
  directed definition accordingly).

**Attributes:** `Protected`.

## References

- D. J. Watts and S. H. Strogatz, *Collective dynamics of 'small-world' networks*, Nature **393** (1998) 440-442.
- N. Chiba and T. Nishizeki, *Arboricity and subgraph listing algorithms*, SIAM J. Comput. **14** (1985) 210-223.
- Source: [`src/graph/gmet_cluster.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_cluster.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

This is the mean of `LocalClusteringCoefficient` over all vertices, where a vertex of degree `d` in `t` triangles scores `t / C(d, 2)`, or 0 when `d < 2`. The answer is exact, an integer or a rational number.

In the triangle-with-pendant example the local values are `1, 1, 1/3, 0`, whose mean is `7/12`. Mixed graphs, with both directed and undirected edges, are left unevaluated.

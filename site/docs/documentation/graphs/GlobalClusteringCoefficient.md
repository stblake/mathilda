# GlobalClusteringCoefficient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GlobalClusteringCoefficient[g] gives 3 x (number of triangles) / (number of connected triples) of g, exactly.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= GlobalClusteringCoefficient[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[1]= 3/5

In[2]:= GlobalClusteringCoefficient[PathGraph[4]]
Out[2]= 0

In[3]:= GlobalClusteringCoefficient[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= 3/4

In[4]:= GlobalClusteringCoefficient[Graph[{1<->2, 2->3}]]
Out[4]= GlobalClusteringCoefficient[Graph[<3 vertices, 2 edges>]]
```

### Applications (3)

Every triple closes, so 1

```mathematica
In[5]:= GlobalClusteringCoefficient[CompleteGraph[4]]
Out[5]= 1
```

No triangles, so 0

```mathematica
In[6]:= GlobalClusteringCoefficient[CycleGraph[5]]
Out[6]= 0
```

The Petersen graph is triangle-free

```mathematica
In[7]:= GlobalClusteringCoefficient[PetersenGraph[]]
Out[7]= 0
```

## Implementation notes

**Algorithm.** `builtin_global_clustering_coefficient` gives
`3 x (triangles) / (connected triples)` exactly — equivalently
`sum t(v) / sum C(d(v), 2)` where `t(v)` is the number of triangles through
vertex `v` and `d(v)` its degree. It shares the triangle machinery of
`GraphTriangleCount`/`LocalClusteringCoefficient`: `tri_compute` lists triangles
over the underlying simple graph using the **degree-ordered orientation** (each
edge points from lower to higher `(degree, index)` rank), which bounds every
oriented out-degree by `O(sqrt m)` and total work by `O(m^1.5)`. For a directed
graph the denominator per vertex is `in(v) out(v) - r(v)` (reciprocal neighbours)
and a "triangle" is a directed 3-cycle. The ratio is returned in lowest terms as
an exact `Rational` (or `0` when the denominator is `0`). A mixed graph is left
unevaluated.

**Data structures.** CSR oriented adjacency (`ooff`/`oadj`, plus a per-edge
arc-flag byte `ofl` when directed), per-thread mark arrays and triangle counters,
and GMP for the final exact ratio. The listing is run on the thread team
(`gmet_parallel_for`).

**Complexity / limits.** `O(m^1.5)` time. Exact (`Integer`/`Rational`);
`EdgeWeight` is ignored. `MeanClusteringCoefficient` is the *mean of the local*
coefficients and is generally a different number from this global one.

- Exact, weights ignored, mixed graphs unevaluated (as in Wolfram).
- Undirected: `3T/Σ C(d, 2)`, with `T` the number of triangles and the sum
  over the vertex degrees `d`.
- Directed (*reverse-engineered*): triangles are directed 3-cycles, with the
  same in/out/reciprocal counting as `LocalClusteringCoefficient`.

**Attributes:** `Protected`.

## References

**See also:** [LocalClusteringCoefficient](../../graphs/LocalClusteringCoefficient/)

- N. Chiba and T. Nishizeki, *Arboricity and subgraph listing algorithms*, SIAM J. Comput. **14** (1985) 210-223.
- Source: [`src/graph/gmet_cluster.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_cluster.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The global (or *transitivity*) coefficient is `3 x (triangles) /
(connected triples)`: the fraction of paths of length two that are closed into
a triangle. It is `1` for a complete graph and `0` for any triangle-free graph.

It is computed exactly as a `Rational`. This is a single graph-wide ratio and
is in general not equal to `MeanClusteringCoefficient`, which averages the
per-vertex local coefficients.

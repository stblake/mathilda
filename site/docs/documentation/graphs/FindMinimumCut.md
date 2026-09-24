# FindMinimumCut

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindMinimumCut[g] gives {value, {part1, part2}}: a partition of the vertices minimizing the total weight of the edges from part1 to part2 (EdgeWeight when present, else 1). For directed graphs part1 is the source side. Nagamochi-Ibaraki for undirected graphs, max flows for directed ones.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= FindMinimumCut[CycleGraph[5]]
Out[1]= {2, {{2, 3, 4, 5}, {1}}}

In[2]:= FindMinimumCut[Graph[{1->2,2->3,3->1,1->3}]]
Out[2]= {1, {{1, 3}, {2}}}

In[3]:= FindMinimumCut[Graph[{1},{}]]
Out[3]= FindMinimumCut[Graph[<1 vertex, 0 edges>]]
```

### Options (2)

```mathematica
In[4]:= FindMinimumCut[Graph[{1,2,3},{UndirectedEdge[1,2],UndirectedEdge[2,3]}, EdgeWeight->{5,2}]]
Out[4]= {2, {{3}, {1, 2}}}

In[5]:= FindMinimumCut[Graph[{1,2,3},{UndirectedEdge[1,2],UndirectedEdge[2,3]}, EdgeWeight->{1/2,2}]]
Out[5]= {0.5, {{2, 3}, {1}}}
```

## Implementation notes

- `Protected`; unevaluated on a non-graph and for fewer than 2 vertices.
- For a graph with directed edges, the edges counted run from `part1` (the
  source side) to `part2`.
- Among equal cuts the shore without `VertexList[g][[1]]` is listed first for
  undirected graphs (Mathematica's tie choice is not reproducible; the value
  always agrees).
- Integer weights give an Integer; Rational or Real weights give a Real (as
  Mathematica); a negative or symbolic weight leaves the call unevaluated.
  Computation is exact in int64 (reals scaled by a common power of two).
- Algorithms: Nagamochi-Ibaraki for undirected global minimum cuts
  (maximum-adjacency orders that contract every edge whose attachment reaches
  the current bound, not one pair per phase as in Stoer-Wagner); `2(n-1)`
  bounded flows (Dinic) for directed global cuts.

**Attributes:** `Protected`.

## References

**See also:** [EdgeWeight](../../graphs/EdgeWeight/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

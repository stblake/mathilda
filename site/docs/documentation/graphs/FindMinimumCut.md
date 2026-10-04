# FindMinimumCut

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindMinimumCut[g] gives {value, {part1, part2}}: a partition of the vertices minimizing the total weight of the edges from part1 to part2 (EdgeWeight when present, else 1). For directed graphs part1 is the source side. Nagamochi-Ibaraki for undirected graphs, max flows for directed ones.`**

## Examples (9)

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

### Applications (4)

One edge separates the path

```mathematica
In[6]:= FindMinimumCut[PathGraph[{1, 2, 3}]]
Out[6]= {1, {{2, 3}, {1}}}
```

A cycle must lose two edges to split

```mathematica
In[7]:= FindMinimumCut[CycleGraph[4]]
Out[7]= {2, {{2, 3, 4}, {1}}}
```

Isolating one vertex of K4 cuts three edges

```mathematica
In[8]:= FindMinimumCut[CompleteGraph[4]]
Out[8]= {3, {{2, 3, 4}, {1}}}
```

Directed: the source side is listed first

```mathematica
In[9]:= FindMinimumCut[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 1 -> 3}]]
Out[9]= {1, {{1, 3}, {2}}}
```

## Implementation notes

**Algorithm.** `builtin_find_minimum_cut` computes a **global** minimum edge cut
of `g` and returns `{value, {part1, part2}}`. Capacities are the `EdgeWeight`
when `g` carries one, else `1`. An all-undirected graph goes to
`gf_ni_mincut`, the Nagamochi-Ibaraki round structure: each round builds a
maximum-adjacency order, lowers the global bound `λ` to the least weighted
degree seen, and contracts every edge whose MA attachment `q(e) ≥ λ` (such edges
carry `λ(u,v) ≥ λ`, so contracting them cannot destroy a lighter cut) together
with the last two vertices of the order. It yields exactly Stoer-Wagner's answer
but typically collapses most of the graph per round. A graph with any directed
edge goes to `gf_directed_mincut`: the minimum over every `v ≠ v₀` of the two
`s-t` max flows `v₀→v` and `v→v₀`, each computed by Dinic's blocking-flow
algorithm (`gf_dinic`) and capped at the running bound.

**Data structures.** Non-integer capacities are scaled by a common power of two
so all arithmetic stays exact on `int64` (a machine real becomes an exact dyadic
integer; a rational goes through its `double`, the precision Mathematica reports).
The flow path uses a CSR residual network (`GfNet`) whose forward and reverse
arcs of one edge are adjacent pairs, with BFS levels and current-arc pointers;
the final cut shore is recovered by residual reachability from the source. The
builtin collects the two shores into `int` arrays and renders them with
`galg_vertex_list` in `VertexList` order.

**Complexity / limits.** Max flow is Dinic (`O(E√E)` on unit capacities, a handful
of phases in practice); the directed global cut runs `O(n)` such flows, the
undirected one is the near-linear-per-round Nagamochi-Ibaraki contraction.
Integer capacities give an exact `Integer` value, rational or real ones a `Real`,
and `Infinity` is an allowed capacity. A graph with fewer than two vertices, or a
negative or symbolic capacity, leaves the call unevaluated; the source side is
listed first for a directed graph, the side without `VertexList[g]⟦1⟧` first for
an undirected one.

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

- H. Nagamochi and T. Ibaraki, *Computing edge-connectivity in multigraphs and capacitated graphs*, SIAM J. Discrete Math. **5** (1992) 54-66.
- M. Stoer and F. Wagner, *A simple min-cut algorithm*, J. ACM **44** (1997) 585-591.
- E. A. Dinic, *Algorithm for solution of a problem of maximum flow in networks with power estimation*, Soviet Math. Dokl. **11** (1970) 1277-1280.
- Source: [`src/graph/galg_flow.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_flow.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`FindMinimumCut[g]` returns `{value, {part1, part2}}`, a global minimum edge cut:
the least total capacity of edges whose removal disconnects `g`, with the two
shores of that cut. Capacities come from `EdgeWeight` when `g` carries one, else
every edge has capacity `1`, so on an unweighted graph the value is the
edge-connectivity.

Undirected graphs use Nagamochi-Ibaraki (the Stoer-Wagner answer, computed by
edge contraction); graphs with a directed edge take the minimum over all vertex
pairs of the two directed max flows. The source side is listed first for a
directed graph, the side without the first vertex first for an undirected one.
Integer capacities give an exact integer value, rational or real ones a real; a
graph on fewer than two vertices, or a negative/symbolic capacity, leaves the
call unevaluated.

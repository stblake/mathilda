# FindPath

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindPath[g, s, t] finds a path from s to t, as {vertex list}, or {} (depth-first, linear time). FindPath[g, s, t, k] finds a path with at most k edges, FindPath[g, s, t, {k}] exactly k, {kmin, kmax} in a range; FindPath[g, s, t, kspec, n] finds at most n paths (n or All), shortest first. Length-bounded searches give up (unevaluated) after 5*10^7 steps.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= FindPath[CycleGraph[5], 1, 3]
Out[1]= {{1, 2, 3}}

In[2]:= FindPath[CompleteGraph[4], 1, 4, 2, All]
Out[2]= {{1, 4}, {1, 3, 4}, {1, 2, 4}}

In[3]:= FindPath[CycleGraph[5], 1, 3, Infinity, All]
Out[3]= {{1, 2, 3}, {1, 5, 4, 3}}

In[4]:= FindPath[Graph[{1->2,3->2}], 1, 3]
Out[4]= {}

In[5]:= FindPath[CycleGraph[5], 2, 2]
Out[5]= {}
```

### Options (1)

```mathematica
In[6]:= FindPath[Graph[{1,2,3},{1<->2,2<->3,3<->1},EdgeWeight->{1,2,3}], 1, 3, 2]
Out[6]= FindPath[Graph[<3 vertices, 3 edges>], 1, 3, 2]
```

## Implementation notes

- `Protected`. A non-graph argument is left unevaluated.
- Plain form: the first path a DFS meets, taking neighbours in `EdgeList` order
  — linear time and the same path as Mathematica.
- `s == t` gives `{}`.
- With a `kspec`, simple paths are enumerated depth-first within the length
  bounds and the first `n` (or `All`) are reported shortest first —
  Mathematica's order.
- Weighted graphs with a `kspec` are left unevaluated (Mathematica measures the
  `kspec` in total weight there).
- Reuses the per-graph incidence-list cache; 1.1–1.6x faster than Mathematica 15
  at `10^5` vertices (`benchmarks/93-graph-ops-editing`).

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

# FindPath

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindPath[g, s, t] finds a path from s to t, as {vertex list}, or {} (depth-first, linear time). FindPath[g, s, t, k] finds a path with at most k edges, FindPath[g, s, t, {k}] exactly k, {kmin, kmax} in a range; FindPath[g, s, t, kspec, n] finds at most n paths (n or All), shortest first. Length-bounded searches give up (unevaluated) after 5*10^7 steps.`**

## Examples (12)

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

### Applications (6)

First path found by depth-first search

```mathematica
In[7]:= FindPath[CycleGraph[6], 1, 4]
Out[7]= {{1, 2, 3, 4}}
```

All paths of exactly three edges

```mathematica
In[8]:= FindPath[CycleGraph[6], 1, 4, {3}, All]
Out[8]= {{1, 6, 5, 4}, {1, 2, 3, 4}}
```

Every simple path, shortest first

```mathematica
In[9]:= FindPath[CompleteGraph[4], 1, 4, Infinity, All]
Out[9]= {{1, 4}, {1, 3, 4}, {1, 2, 4}, {1, 3, 2, 4}, {1, 2, 3, 4}}
```

No directed path against the arrows

```mathematica
In[10]:= FindPath[Graph[{1 -> 2, 2 -> 3}], 3, 1]
Out[10]= {}
```

Disconnected components give no path

```mathematica
In[11]:= FindPath[Graph[{1 <-> 2, 3 <-> 4}], 1, 4]
Out[11]= {}
```

Start equal to end gives an empty result

```mathematica
In[12]:= FindPath[CycleGraph[6], 1, 1]
Out[12]= {}
```

## Implementation notes

**Algorithm.** `builtin_find_path` takes `FindPath[g, s, t]`, `FindPath[g, s, t, kspec]` and `FindPath[g, s, t, kspec, n]`, with `kspec` as for `FindCycle` (length counts edges) and `n` a positive integer or `All`. The basic form is an iterative DFS from `s` (`dfs_path`) with visited marks and neighbours taken in `EdgeList` order, returning `{path}` as a vertex list, or `{}` when `t` is unreachable. `s == t` gives `{}`. Bounded or multi-path forms backtrack over simple paths with an on-path mark, record every path that reaches `t` within the length range, and then order them shortest first (ties by discovery order).

**Data structures.** A cached forward-arc CSR built from the edge-index views (undirected edges contribute both arcs); `int` arrays for the DFS stack of (vertex, cursor) pairs, the current path and the on-path marks. The result is a list of vertex lists.

**Complexity / limits.** `O(V + E)` for the single-path form, and a repeated query on one graph reuses the cached CSR. Path enumeration is exponential in the worst case and is bounded by a `TimeConstrained`-aware step budget, after which the call stays unevaluated. A vertex not in the graph, or a weighted graph with a `kspec`, is also left unevaluated.

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

- Source: [`src/graph/gops_cycles.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_cycles.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

A path is a list of vertices and the result is a list of paths, `{}` if none exists. Length counts edges, with the same `kspec` forms as `FindCycle`, and a fifth argument bounds the number of paths (or `All`); multiple paths are listed shortest first.

The single-path search is a linear-time depth-first search with neighbours taken in `EdgeList` order. Enumeration is bounded by a step budget, past which the call stays unevaluated.

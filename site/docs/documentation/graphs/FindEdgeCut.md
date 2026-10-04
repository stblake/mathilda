# FindEdgeCut

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindEdgeCut[g] gives a minimum set of edges whose removal disconnects g (strongly, for directed graphs); FindEdgeCut[g, s, t] gives a minimum s-t edge cut. Uses EdgeWeight as capacity when present. Edges are returned in EdgeList order.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= FindEdgeCut[CycleGraph[5]]
Out[1]= {1 <-> 2, 5 <-> 1}

In[2]:= FindEdgeCut[CycleGraph[6], 1, 4]
Out[2]= {1 <-> 2, 6 <-> 1}

In[3]:= FindEdgeCut[Graph[{1->2,2->3,3->1}]]
Out[3]= {1 -> 2}

In[4]:= FindEdgeCut[x]
Out[4]= FindEdgeCut[x]
```

### Options (1)

```mathematica
In[5]:= FindEdgeCut[Graph[{1,2,3},{UndirectedEdge[1,2],UndirectedEdge[2,3],UndirectedEdge[3,1]}, EdgeWeight->{2,3,4}]]
Out[5]= {1 <-> 2, 2 <-> 3}
```

### Applications (3)

A single bridge is enough

```mathematica
In[6]:= FindEdgeCut[PathGraph[{1, 2, 3}]]
Out[6]= {1 <-> 2}
```

Two edges, here both incident to one vertex

```mathematica
In[7]:= FindEdgeCut[CycleGraph[4]]
Out[7]= {1 <-> 2, 4 <-> 1}
```

A minimum 1-2 edge cut

```mathematica
In[8]:= FindEdgeCut[CompleteGraph[4], 1, 2]
Out[8]= {1 <-> 2, 1 <-> 3, 1 <-> 4}
```

## Implementation notes

**Algorithm.** `builtin_find_edge_cut` returns the *edges* of a minimum cut (in
`EdgeList` order), rather than its weight. `FindEdgeCut[g]` finds a global
minimum cut — Nagamochi-Ibaraki for an undirected graph, the min over `s-t`
flows for a directed one — then reports the edges crossing the shore boundary
(`gf_cut_edges`: an undirected edge whenever its endpoints lie on opposite
sides, a directed one only when it runs side-1 to side-0). `FindEdgeCut[g, s, t]`
runs a single Dinic max flow from `s` to `t` and takes the cut closest to `s`
(the set of vertices still residual-reachable from `s`). `EdgeWeight` is used as
the capacity when present, else `1`.

**Data structures.** The same `GfCap` scaled-capacity parse and `GfNet` CSR
residual network as `EdgeConnectivity` and `FindMaximumFlow`; the shore is a
`char side[]` filled by a residual-reachability BFS (`gf_reach_from`), from which
the crossing edges are collected over the memo's integer endpoint arrays.

**Complexity / limits.** Dominated by the flow/contraction phase. The `s-t` form
returns `{}` when `s` and `t` are already separated by an empty cut; `FindEdgeCut`
of a graph with one vertex is `{}`. A negative or symbolic capacity leaves the
call unevaluated.

- `Protected`; unevaluated on a non-graph.
- Cuts are weighted by `EdgeWeight`; the edges are returned in `EdgeList` order.
- The `s`-`t` cut is the one closest to `s`, as in Mathematica.
- Directed graphs use strong connectivity.
- Uses the same exact int64 flow machinery as `FindMaximumFlow` /
  `FindMinimumCut` (Dinic; Nagamochi-Ibaraki for undirected global cuts).

**Attributes:** `Protected`.

## References

**See also:** [EdgeWeight](../../graphs/EdgeWeight/), [EdgeList](../../graphs/EdgeList/), [FindMaximumFlow](../../graphs/FindMaximumFlow/), [FindMinimumCut](../../graphs/FindMinimumCut/)

- H. Nagamochi and T. Ibaraki, *Computing edge-connectivity in multigraphs and capacitated graphs*, SIAM J. Discrete Math. **5** (1992) 54-66.
- E. A. Dinic, *Algorithm for solution of a problem of maximum flow in a network with power estimation*, Soviet Math. Dokl. **11** (1970) 1277-1280.
- Source: [`src/graph/galg_flow.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_flow.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`FindEdgeCut[g]` returns the edges of a global minimum cut; `FindEdgeCut[g, s,
t]` the edges of a minimum cut separating `s` from `t`. The number of edges (or
their total `EdgeWeight`) equals `EdgeConnectivity` of the same arguments.

Edges are returned in `EdgeList` order. The `s-t` cut is the one closest to
`s`; adjacent `s` and `t` on an unweighted graph can only be separated by
removing every edge between them.

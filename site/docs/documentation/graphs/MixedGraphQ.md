# MixedGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MixedGraphQ[g] gives True if g has both directed and undirected edges.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= {MixedGraphQ[Graph[{1->2, 2<->3}]], MixedGraphQ[Graph[{1->2}]], MixedGraphQ[CycleGraph[3]], MixedGraphQ[x]}
Out[1]= {True, False, False, False}
```

### Applications (3)

One directed and one undirected edge

```mathematica
In[2]:= MixedGraphQ[Graph[{1, 2, 3}, {1 -> 2, 2 <-> 3}]]
Out[2]= True
```

All undirected

```mathematica
In[3]:= MixedGraphQ[CycleGraph[3]]
Out[3]= False
```

All directed

```mathematica
In[4]:= MixedGraphQ[Graph[{1, 2}, {1 -> 2}]]
Out[4]= False
```

## Implementation notes

**Algorithm.** `builtin_mixed_graph_q` is `True` exactly when `g` has **both** a
directed and an undirected edge. It reads the directed-edge count from the memo
(`graph_directed_edge_count`) and compares it with the total edge count: the
graph is mixed iff the count is strictly between `0` and the total. A graph with
no directed edges (all undirected or edgeless) or with every edge directed is not
mixed.

**Data structures.** None of its own — two integers off the validated-graph memo
entry. The result is a fresh `True`/`False` symbol.

**Complexity / limits.** `O(1)` on a memoized graph. A non-graph argument gives
`False`, never unevaluated — as every `*Q` predicate does. Several transforms
(`UndirectedGraph`, `DirectedGraph`, `LineGraph`, `FindSpanningTree`) decline on
mixed graphs, so this is the predicate that tells them apart.

- `Protected`. A non-graph argument gives `False`.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_preds.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_preds.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

A graph is *mixed* when it has at least one directed edge and at least one
undirected edge. A purely directed, purely undirected, or edgeless graph is not
mixed.

Several operations — `UndirectedGraph`, `DirectedGraph`, `LineGraph`,
`FindSpanningTree`, the clustering coefficients — are left unevaluated on mixed
graphs, so `MixedGraphQ` is the guard for them.

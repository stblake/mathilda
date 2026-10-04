# ReverseGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ReverseGraph[g] reverses every directed edge of g; undirected edges, the edge order and the weights are kept.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= EdgeList[ReverseGraph[Graph[{1->2,2->3,3<->4}]]]
Out[1]= {2 -> 1, 3 -> 2, 3 <-> 4}

In[2]:= EdgeList[ReverseGraph[CycleGraph[3]]]
Out[2]= {1 <-> 2, 2 <-> 3, 3 <-> 1}
```

### Options (1)

```mathematica
In[3]:= InputForm[ReverseGraph[Graph[{1,2,3},{1->2,3->2},EdgeWeight->{4,5}]]]
Out[3]= Graph[{1, 2, 3}, {2 -> 1, 2 -> 3}, EdgeWeight -> {4, 5}]
```

### Applications (2)

Each arc flips direction

```mathematica
In[4]:= EdgeList[ReverseGraph[Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}]]]
Out[4]= {2 -> 1, 3 -> 2}
```

Undirected edges are unchanged

```mathematica
In[5]:= EdgeList[ReverseGraph[CycleGraph[3]]]
Out[5]= {1 <-> 2, 2 <-> 3, 3 <-> 1}
```

## Implementation notes

**Algorithm.** `builtin_reverse_graph` reverses every directed edge of `g` while
leaving undirected edges, the edge order, and the weights untouched. When `g` has
no directed edges it returns a copy of the input unchanged (a fast exit). For a
graph with directed edges it walks the edge list and, for each `DirectedEdge[a,
b]`, emits `DirectedEdge[b, a]` with its endpoints swapped; an undirected edge is
copied verbatim. The `k`-th output edge corresponds to the `k`-th input edge, so
`EdgeList` order is preserved and the `EdgeWeight`/`EdgeCapacity` lists stay
aligned.

**Data structures.** A `GopsView` of the input and fresh vertex, edge, weight,
and integer-endpoint arrays; reversed edges share a single interned
`DirectedEdge` head (`GopsHeads`). `gops_graph_new` seeds the result into the
graph memo.

**Complexity / limits.** `O(V + E)`, a single linear pass. Weights are carried
over unchanged, so a weighted reversed graph keeps its weights. A non-graph
argument leaves the call unevaluated.

- `Protected`. A non-graph argument is left unevaluated.
- Undirected edges are kept as they are; edge order and weights are kept.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_transform.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_transform.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`ReverseGraph[g]` reverses the direction of every directed edge; undirected
edges, the edge order, and any `EdgeWeight` are kept. An undirected graph is
returned unchanged.

Reversing twice recovers the original graph.

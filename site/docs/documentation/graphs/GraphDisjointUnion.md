# GraphDisjointUnion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphDisjointUnion[g1, g2, ...] gives the disjoint union of the gi, with vertices relabelled 1, 2, ..., n (g1's first, in VertexList order). Weights are dropped.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= EdgeList[GraphDisjointUnion[Graph[{a<->b}], Graph[{a<->b, b<->c}]]]
Out[1]= {1 <-> 2, 3 <-> 4, 4 <-> 5}

In[2]:= VertexList[GraphDisjointUnion[CycleGraph[3], PathGraph[{x,y}]]]
Out[2]= {1, 2, 3, 4, 5}

In[3]:= InputForm[GraphDisjointUnion[Graph[{a<->b}]]]
Out[3]= Graph[{a, b}, {a <-> b}]
```

### Applications (3)

Vertices never merge: 3 + 4

```mathematica
In[4]:= VertexCount[GraphDisjointUnion[CycleGraph[3], CycleGraph[4]]]
Out[4]= 7
```

3 + 3 edges, kept apart

```mathematica
In[5]:= EdgeCount[GraphDisjointUnion[CompleteGraph[3], CompleteGraph[3]]]
Out[5]= 6
```

Relabelled 1..n

```mathematica
In[6]:= VertexList[GraphDisjointUnion[CycleGraph[3], PathGraph[{1, 2}]]]
Out[6]= {1, 2, 3, 4, 5}
```

## Implementation notes

**Algorithm.** `builtin_graph_disjoint_union` places the input graphs side by
side on fresh integer vertices `1..n`: graph `g1`'s vertices first (in
`VertexList` order), then `g2`'s, and so on, with each graph's edges translated
by its block offset. No vertices are identified even when the inputs share labels
— this is exactly what distinguishes it from `GraphUnion`, which merges equal
vertices. `GraphDisjointUnion[g]` is `g`. Every edge keeps its direction; weights
are dropped.

**Data structures.** It first sums the vertex and edge counts to size the output
arrays, fills the relabelled vertex list `{1, ..., n}`, then opens a `GopsView`
on each input and copies its edges with endpoints shifted by a running `off`.
Edge nodes are built through a shared interned-head cache (`GopsHeads`) and the
result goes through `gops_graph_new`, which seeds the graph memo.

**Complexity / limits.** `O(V + E)` over all inputs — a single linear pass, no
hashing (unlike the vertex-union set operations). Each argument must be a valid
graph, else the call is unevaluated.

- `Protected`. A non-graph argument is left unevaluated.
- Vertices are relabelled `1..n`, `g1`'s first; edges are translated in order.
- `GraphDisjointUnion[g]` is `g` (not relabelled).
- Weights are dropped, as in Mathematica. Shares the set-operation machinery
  described under `GraphUnion`.

**Attributes:** `Protected`.

## References

**See also:** [GraphUnion](../../graphs/GraphUnion/)

- Source: [`src/graph/gops_setops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_setops.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
- Tests: [`tests/test_graphplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphplot.c)

## Notes & additional examples

### Notes

`GraphDisjointUnion` lays the graphs out on disjoint vertex sets even when they
share vertex names, relabelling everything to `1, ..., n` with the first graph's
vertices first. The result's vertex and edge counts are the plain sums of the
inputs'.

Contrast `GraphUnion`, which identifies equal vertices across its arguments.
Edge weights are dropped.

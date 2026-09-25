# GraphDisjointUnion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphDisjointUnion[g1, g2, ...] gives the disjoint union of the gi, with vertices relabelled 1, 2, ..., n (g1's first, in VertexList order). Weights are dropped.`**

## Examples (3)

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

## Implementation notes

- `Protected`. A non-graph argument is left unevaluated.
- Vertices are relabelled `1..n`, `g1`'s first; edges are translated in order.
- `GraphDisjointUnion[g]` is `g` (not relabelled).
- Weights are dropped, as in Mathematica. Shares the set-operation machinery
  described under `GraphUnion`.

**Attributes:** `Protected`.

## References

**See also:** [GraphUnion](../../graphs/GraphUnion/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

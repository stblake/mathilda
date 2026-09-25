# VertexCoverQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexCoverQ[g, vlist] gives True if every edge of g has an endpoint in vlist.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= IndependentVertexSetQ[CycleGraph[6], {1, 3, 5}]
Out[1]= True

In[2]:= IndependentVertexSetQ[CycleGraph[6], {1, 2}]
Out[2]= False

In[3]:= IndependentVertexSetQ[CycleGraph[6], {1, 1, 3}]
Out[3]= True

In[4]:= IndependentVertexSetQ[CycleGraph[6], {1, 7}]
Out[4]= False

In[5]:= VertexCoverQ[CycleGraph[4], {1, 3}]
Out[5]= True

In[6]:= VertexCoverQ[x, {1}]
Out[6]= False
```

## Implementation notes

- `Protected` membership predicates; give `False` (never stay unevaluated) when
  `g` is not a graph.
- An element that is not a vertex of `g` gives `False`.
- Repeated vertices are allowed in `vs`.

**Attributes:** `Protected`.

## References

**See also:** [IndependentVertexSetQ](../../graphs/IndependentVertexSetQ/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

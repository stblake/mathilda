# FindGraphIsomorphism

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindGraphIsomorphism[g1, g2] gives {assoc}, an isomorphism from g1 to g2 as an association of vertices, or {} if the graphs are not isomorphic. FindGraphIsomorphism[g1, g2, n] / [g1, g2, All] gives up to n / all isomorphisms.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= FindGraphIsomorphism[CycleGraph[4], Graph[{a,b,c,d},{UndirectedEdge[a,c],UndirectedEdge[c,b],UndirectedEdge[b,d],UndirectedEdge[d,a]}]]
Out[1]= {<|1 -> a, 2 -> c, 3 -> b, 4 -> d|>}

In[2]:= FindGraphIsomorphism[CycleGraph[4], CycleGraph[4], 2]
Out[2]= {<|1 -> 1, 2 -> 2, 3 -> 3, 4 -> 4|>, <|1 -> 1, 2 -> 4, 3 -> 3, 4 -> 2|>}

In[3]:= Length[FindGraphIsomorphism[CycleGraph[4], CycleGraph[4], All]]
Out[3]= 8

In[4]:= FindGraphIsomorphism[CycleGraph[4], StarGraph[4]]
Out[4]= {}

In[5]:= FindGraphIsomorphism[Graph[{},{}], Graph[{},{}]]
Out[5]= {<||>}
```

## Implementation notes

- `Protected`; unevaluated on non-graph arguments.
- The Association runs over `VertexList[g1]` in order.
- The order of the list of isomorphisms is the engine's search order;
  Mathematica's differs.
- Two empty graphs give `{<||>}` (Mathematica gives `{}`, which contradicts its
  own `IsomorphicGraphQ` answer `True`).
- Uses the `IsomorphicGraphQ` individualization-refinement engine; every map
  returned is verified edge by edge. Directed and mixed graphs are supported;
  edge weights are ignored.

**Attributes:** `Protected`.

## References

**See also:** [IsomorphicGraphQ](../../graphs/IsomorphicGraphQ/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

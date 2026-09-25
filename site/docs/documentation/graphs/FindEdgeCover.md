# FindEdgeCover

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindEdgeCover[g] gives a minimum edge cover of g: a smallest set of edges touching every vertex. Gives {} when g has an isolated vertex (no edge cover exists).`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= FindEdgeCover[PathGraph[{1,2,3,4,5}]]
Out[1]= {1 <-> 2, 2 <-> 3, 4 <-> 5}

In[2]:= FindEdgeCover[StarGraph[4]]
Out[2]= {1 <-> 2, 1 <-> 3, 1 <-> 4}

In[3]:= FindEdgeCover[Graph[{1,2,3},{UndirectedEdge[1,2]}]]
Out[3]= {}
```

## Implementation notes

- `Protected`; unevaluated on a non-graph.
- Computed as a maximum matching (`FindIndependentEdgeSet`) plus one edge per
  exposed vertex.
- `{}` when `g` has an isolated vertex, as in Mathematica.
- Weighted graphs: optimizes cardinality, as the Wolfram documentation states.
  (The differential test found Mathematica returning non-minimum edge covers on
  weighted graphs.)

**Attributes:** `Protected`.

## References

**See also:** [FindIndependentEdgeSet](../../graphs/FindIndependentEdgeSet/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

# EdgeCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeCount[g] gives the number of edges in the graph g.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= VertexCount[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[1]= 4

In[2]:= EdgeCount[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[2]= 4

In[3]:= VertexCount[Graph[{}, {}]]
Out[3]= 0

In[4]:= EdgeCount[Graph[{1,2,3},{}]]
Out[4]= 0

In[5]:= EdgeCount[x]
Out[5]= EdgeCount[x]
```

## Implementation notes

- `Protected`. Cardinalities read from the canonical form; unevaluated on a
  non-graph (see `VertexList`).

**Attributes:** `Protected`.

## References

**See also:** [VertexCount](../../graphs/VertexCount/), [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

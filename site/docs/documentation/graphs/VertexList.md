# VertexList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexList[g] gives the list of vertices of the graph g.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= VertexList[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[1]= {1, 2, 3, 4}

In[2]:= VertexList[Graph[{c->a, a->b}]]
Out[2]= {c, a, b}

In[3]:= VertexList[Graph[{}, {}]]
Out[3]= {}

In[4]:= VertexList[5]
Out[4]= VertexList[5]
```

## Implementation notes

- `Protected`. The query/representation heads (`VertexList`, `EdgeList`,
  `VertexCount`, `EdgeCount`, `AdjacencyList`, `VertexDegree`,
  `VertexInDegree`, `VertexOutDegree`, `EdgeWeight`) are all thin readers over
  the canonical form and return unevaluated on a non-graph argument.
- Canonical order is the explicit vertex list, or first-appearance order when
  `Graph[e]` derived the vertices from the edges.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/), [VertexCount](../../graphs/VertexCount/), [EdgeCount](../../graphs/EdgeCount/), [AdjacencyList](../../graphs/AdjacencyList/), [VertexDegree](../../graphs/VertexDegree/), [VertexInDegree](../../graphs/VertexInDegree/), [VertexOutDegree](../../graphs/VertexOutDegree/), [EdgeWeight](../../graphs/EdgeWeight/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

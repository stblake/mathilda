# EdgeList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeList[g] gives the list of edges of the graph g.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[Graph[{1,2,3},{1->2, 2<->3}]]
Out[1]= {1 -> 2, 2 <-> 3}

In[2]:= EdgeList[CycleGraph[4]]
Out[2]= {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}

In[3]:= EdgeList[Graph[{1,2},{}]]
Out[3]= {}

In[4]:= EdgeList[5]
Out[4]= EdgeList[5]
```

## Implementation notes

- `Protected`. A thin reader over the canonical form; unevaluated on a non-graph
  (see `VertexList`).
- Edges print in operator form (`1 -> 2`, `2 <-> 3`), but are
  `DirectedEdge`/`UndirectedEdge` internally.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

# PathGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PathGraph[n] gives the path on n vertices; PathGraph[{v1,...}] the path over the given vertices.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= VertexDegree[PathGraph[5]]
Out[1]= {1, 2, 2, 2, 1}

In[2]:= EdgeList[PathGraph[4]]
Out[2]= {1 <-> 2, 2 <-> 3, 3 <-> 4}

In[3]:= EdgeList[PathGraph[{a,b,c}]]
Out[3]= {a <-> b, b <-> c}

In[4]:= PathGraph[1]
Out[4]= Graph[<1 vertex, 0 edges>]

In[5]:= PathGraph[x]
Out[5]= PathGraph[x]
```

## Implementation notes

- `Protected`. Undirected edges, built through the `Graph` constructor (see
  `CycleGraph`). `PathGraph[1]` is a single vertex; a symbolic argument is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/), [CycleGraph](../../graphs/CycleGraph/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

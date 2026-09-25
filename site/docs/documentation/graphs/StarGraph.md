# StarGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StarGraph[n] gives the star on n vertices: the hub 1 joined to each of the n-1 leaves 2..n.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[StarGraph[4]]
Out[1]= {1 <-> 2, 1 <-> 3, 1 <-> 4}

In[2]:= VertexDegree[StarGraph[5]]
Out[2]= {4, 1, 1, 1, 1}

In[3]:= StarGraph[1]
Out[3]= Graph[<1 vertex, 0 edges>]

In[4]:= StarGraph[x]
Out[4]= StarGraph[x]
```

## Implementation notes

- `Protected`. Undirected edges, built through the `Graph` constructor (see
  `CycleGraph`). `StarGraph[1]` is a single vertex; a symbolic argument is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/), [CycleGraph](../../graphs/CycleGraph/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

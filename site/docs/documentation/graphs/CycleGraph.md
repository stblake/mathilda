# CycleGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CycleGraph[n] gives the cycle graph on n vertices.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= EdgeList[CycleGraph[4]]
Out[1]= {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}

In[2]:= VertexCount[CycleGraph[10]]
Out[2]= 10

In[3]:= EdgeList[CycleGraph[2]]
Out[3]= {1 <-> 2}

In[4]:= EdgeCount[CompleteGraph[5]]
Out[4]= 10

In[5]:= CycleGraph[x]
Out[5]= CycleGraph[x]
```

## Implementation notes

- `Protected`. Like all the graph generators (`CompleteGraph`, `CycleGraph`,
  `PathGraph`, `StarGraph`, `RandomGraph`), it builds a canonical graph with
  vertices `1..n` and undirected edges via the `Graph` constructor path.
- Related generator `CompleteGraph[n]` — `K_n`, with all `n(n-1)/2` edges.
- Small cases: `CycleGraph[2]` is the single edge `1 <-> 2` (no parallel edges),
  `CycleGraph[1]` is one isolated vertex and `CycleGraph[0]` the null graph. A
  symbolic `n` is left unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [CompleteGraph](../../graphs/CompleteGraph/), [PathGraph](../../graphs/PathGraph/), [StarGraph](../../graphs/StarGraph/), [RandomGraph](../../graphs/RandomGraph/), [Graph](../../graphs/Graph/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

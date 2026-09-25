# CompleteGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CompleteGraph[n] gives the complete graph on n vertices. CompleteGraph[{n1, n2, ...}] gives the complete multipartite graph with parts of sizes n1, n2, ....`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= CompleteGraph[5]
Out[1]= Graph[<5 vertices, 10 edges>]

In[2]:= CompleteGraph[{2, 3}]
Out[2]= Graph[<5 vertices, 6 edges>]

In[3]:= EdgeList[CompleteGraph[{2, 3}]]
Out[3]= {1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 4, 2 <-> 5}

In[4]:= EdgeList[CompleteGraph[{1, 1, 2}]]
Out[4]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 3, 2 <-> 4}

In[5]:= CompleteGraph[{4}]
Out[5]= Graph[<4 vertices, 6 edges>]

In[6]:= CompleteGraph[x]
Out[6]= CompleteGraph[x]
```

## Implementation notes

- All graph families of this module are undirected on `1..N`; the edge list is
  the sorted list of pairs `{i, j}`, `i < j` — identical to Wolfram's
  `EdgeList` for each family.
- `CompleteGraph[{n}]` is `K_n`.
- `CompleteGraph` is re-registered by a wrapper that delegates its
  pre-existing form (`CompleteGraph[n]`) to the original builtin.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: a family with more than 10^8 vertices or 5×10^7 edges is
  left unevaluated rather than allocating gigabytes.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

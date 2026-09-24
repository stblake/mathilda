# TuranGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TuranGraph[n, k] gives the Turan graph: the complete k-partite graph on n vertices with parts as equal as possible.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= TuranGraph[5, 2]
Out[1]= Graph[<5 vertices, 6 edges>]

In[2]:= EdgeList[TuranGraph[5, 2]]
Out[2]= {1 <-> 4, 1 <-> 5, 2 <-> 4, 2 <-> 5, 3 <-> 4, 3 <-> 5}

In[3]:= TuranGraph[6, 3]
Out[3]= Graph[<6 vertices, 12 edges>]

In[4]:= EdgeList[TuranGraph[4, 3]]
Out[4]= {1 <-> 3, 1 <-> 4, 2 <-> 3, 2 <-> 4, 3 <-> 4}
```

## Implementation notes

- Undirected on `1..n`, larger parts first; the edge list is the sorted list
  of pairs `{i, j}`, `i < j`, identical to Wolfram's `EdgeList`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

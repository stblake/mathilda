# WheelGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WheelGraph[n] gives the wheel graph with n vertices: vertex 1 joined to every vertex of the cycle 2, ..., n.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= WheelGraph[5]
Out[1]= Graph[<5 vertices, 8 edges>]

In[2]:= EdgeList[WheelGraph[5]]
Out[2]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 5, 3 <-> 4, 4 <-> 5}

In[3]:= WheelGraph[1]
Out[3]= Graph[<1 vertex, 0 edges>]

In[4]:= WheelGraph[3]
Out[4]= WheelGraph[3]
```

### Options (1)

```mathematica
In[5]:= WheelGraph[5, DirectedEdges -> True]
Out[5]= WheelGraph[5, DirectedEdges -> True]
```

## Implementation notes

- Undirected on `1..n`; the edge list is the sorted list of pairs `{i, j}`,
  `i < j`, identical to Wolfram's `EdgeList`.
- Accepts `n = 1` or `n >= 4`; `n = 2` and `3` would be multigraphs and are
  left unevaluated.
- Options (`DirectedEdges`, layout options) are not supported: a call with
  options is left unevaluated.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

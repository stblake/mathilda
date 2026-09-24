# HararyGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HararyGraph[k, n] gives the Harary graph: a k-connected graph on n vertices with the minimum number of edges (k >= 2, n > k).`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= HararyGraph[2, 5]
Out[1]= Graph[<5 vertices, 5 edges>]

In[2]:= EdgeList[HararyGraph[3, 6]]
Out[2]= {1 <-> 2, 1 <-> 4, 1 <-> 6, 2 <-> 3, 2 <-> 5, 3 <-> 4, 3 <-> 6, 4 <-> 5, 5 <-> 6}

In[3]:= EdgeList[HararyGraph[3, 5]]
Out[3]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 5, 3 <-> 4, 4 <-> 5}

In[4]:= HararyGraph[1, 5]
Out[4]= HararyGraph[1, 5]

In[5]:= HararyGraph[4, 4]
Out[5]= HararyGraph[4, 4]
```

## Implementation notes

- Requires `k >= 2` and `n > k`; otherwise unevaluated.
- Undirected on `1..n`; the edge list is the sorted list of pairs `{i, j}`,
  `i < j`, identical to Wolfram's `EdgeList`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

# CirculantGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CirculantGraph[n, j] gives the circulant graph on n vertices with i joined to i+j and i-j (mod n); CirculantGraph[n, {j1, j2, ...}] uses every jump ji.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= CirculantGraph[5, 1]
Out[1]= Graph[<5 vertices, 5 edges>]

In[2]:= EdgeList[CirculantGraph[6, 2]]
Out[2]= {1 <-> 3, 1 <-> 5, 2 <-> 4, 2 <-> 6, 3 <-> 5, 4 <-> 6}

In[3]:= EdgeList[CirculantGraph[6, {1, 3}]]
Out[3]= {1 <-> 2, 1 <-> 4, 1 <-> 6, 2 <-> 3, 2 <-> 5, 3 <-> 4, 3 <-> 6, 4 <-> 5, 5 <-> 6}

In[4]:= VertexDegree[CirculantGraph[8, {1, 2}]]
Out[4]= {4, 4, 4, 4, 4, 4, 4, 4}
```

## Implementation notes

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

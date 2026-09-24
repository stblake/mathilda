# PetersenGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PetersenGraph[] gives the Petersen graph; PetersenGraph[n, k] the generalized Petersen graph: inner vertices 1..n with i joined to i+k (mod n), outer cycle n+1..2n, and spokes i to n+i.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= PetersenGraph[]
Out[1]= Graph[<10 vertices, 15 edges>]

In[2]:= VertexDegree[PetersenGraph[]]
Out[2]= {3, 3, 3, 3, 3, 3, 3, 3, 3, 3}

In[3]:= GraphDiameter[PetersenGraph[]]
Out[3]= 2

In[4]:= EdgeList[PetersenGraph[4, 1]]
Out[4]= {1 <-> 2, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 6, 3 <-> 4, 3 <-> 7, 4 <-> 8, 5 <-> 6, 5 <-> 8, 6 <-> 7, 7 <-> 8}
```

## Implementation notes

- Undirected on `1..2n`: the inner star is `1..n`, the outer cycle `n+1..2n`;
  the edge list is the sorted list of pairs `{i, j}`, `i < j`, identical to
  Wolfram's `EdgeList`.
- `PetersenGraph[]` is `PetersenGraph[5, 2]`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

# GridGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GridGraph[{m, n}] gives the m x n grid graph; GridGraph[{n1, ..., nk}] the k-dimensional grid. The first coordinate varies fastest in the vertex numbering.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= GridGraph[{3, 2}]
Out[1]= Graph[<6 vertices, 7 edges>]

In[2]:= EdgeList[GridGraph[{3, 2}]]
Out[2]= {1 <-> 2, 1 <-> 4, 2 <-> 3, 2 <-> 5, 3 <-> 6, 4 <-> 5, 5 <-> 6}

In[3]:= GridGraph[{2, 2, 2}]
Out[3]= Graph[<8 vertices, 12 edges>]

In[4]:= GridGraph[{5}]
Out[4]= Graph[<5 vertices, 4 edges>]

In[5]:= GridGraph[Table[2, {26}]]
Out[5]= GridGraph[{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}]
```

## Implementation notes

- Undirected on `1..N`, vertices numbered with the first coordinate varying
  fastest; the edge list is the sorted list of pairs `{i, j}`, `i < j`,
  identical to Wolfram's `EdgeList`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limits (shared by every graph family of this module): a family
  with more than 10^8 vertices or 5×10^7 edges is left unevaluated rather than
  allocating gigabytes (e.g. `GridGraph[Table[2, {26}]]`, ~8.7×10^8 edges).

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

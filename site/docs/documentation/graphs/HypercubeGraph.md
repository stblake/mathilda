# HypercubeGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypercubeGraph[n] gives the n-dimensional hypercube graph on 2^n vertices.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= HypercubeGraph[3]
Out[1]= Graph[<8 vertices, 12 edges>]

In[2]:= EdgeList[HypercubeGraph[2]]
Out[2]= {1 <-> 2, 1 <-> 3, 2 <-> 4, 3 <-> 4}

In[3]:= VertexDegree[HypercubeGraph[4]]
Out[3]= {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}

In[4]:= HypercubeGraph[0]
Out[4]= Graph[<1 vertex, 0 edges>]
```

## Implementation notes

- Undirected on `1..2^n`; the edge list is the sorted list of pairs `{i, j}`,
  `i < j`, identical to Wolfram's `EdgeList`.
- `HypercubeGraph[0]` is a single vertex.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

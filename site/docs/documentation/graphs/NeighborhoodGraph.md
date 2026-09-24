# NeighborhoodGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NeighborhoodGraph[g, v] gives the subgraph of g induced by v and its neighbours; NeighborhoodGraph[g, v, k] by the vertices within distance k of v (k a non-negative integer or Infinity), edge direction ignored. NeighborhoodGraph[g, {v1, ...}, k] uses several centres. Vertices: the centres, then each centre's new vertices in VertexList order.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= VertexList[NeighborhoodGraph[PathGraph[Range[6]], {1,6}]]
Out[1]= {1, 6, 2, 5}

In[2]:= VertexList[NeighborhoodGraph[PathGraph[Range[6]], 3, 2]]
Out[2]= {3, 1, 2, 4, 5}

In[3]:= EdgeList[NeighborhoodGraph[Graph[{1->2,3->1,2->4}], 1]]
Out[3]= {1 -> 2, 3 -> 1}

In[4]:= VertexList[NeighborhoodGraph[PathGraph[Range[6]], 2, Infinity]]
Out[4]= {2, 1, 3, 4, 5, 6}

In[5]:= VertexList[NeighborhoodGraph[PathGraph[Range[6]], 9]]
Out[5]= {}
```

## Implementation notes

- `Protected`. A non-graph first argument is left unevaluated.
- Distance ignores edge direction.
- Vertices: the centres, then each centre's new vertices in `VertexList` order;
  edges are ordered as for `Subgraph`.
- Non-vertex centres are ignored.
- `k = Infinity` is accepted (Mathematica leaves it unevaluated).

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/), [Subgraph](../../graphs/Subgraph/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

# VertexDelete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexDelete[g, v] removes the vertex v and its incident edges from g; VertexDelete[g, {v1, ...}] removes several (each must be a vertex of g, else the call stays unevaluated); VertexDelete[g, patt] removes every vertex matching patt. Vertex and edge order and edge weights are kept.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= EdgeList[VertexDelete[CycleGraph[5], {1, 3}]]
Out[1]= {4 <-> 5}

In[2]:= VertexList[VertexDelete[PathGraph[Range[6]], _?EvenQ]]
Out[2]= {1, 3, 5}

In[3]:= VertexDelete[CycleGraph[3], 7]
Out[3]= VertexDelete[Graph[<3 vertices, 3 edges>], 7]
```

### Options (1)

```mathematica
In[4]:= InputForm[VertexDelete[Graph[{1,2,3,4},{1<->2,2<->3,3<->4},EdgeWeight->{5,6,7}], 2]]
Out[4]= Graph[{1, 3, 4}, {3 <-> 4}, EdgeWeight -> {7}]
```

## Implementation notes

- `Protected`. A non-graph first argument is left unevaluated.
- Every listed vertex must exist; otherwise the call is left unevaluated.
- Orders and weights of the surviving vertices and edges are kept.
- `O(V + E)` integer pass over the memoized endpoint arrays plus one hash per
  argument item; the result is seeded into the graph memo (see `VertexAdd`).

**Attributes:** `Protected`.

## References

**See also:** [VertexAdd](../../graphs/VertexAdd/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

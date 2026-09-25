# UndirectedGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UndirectedGraph[g] gives the undirected graph underlying g: u->v and v->u become the single edge u<->v, whose weight is the sum of theirs. Edges are oriented and ordered by VertexList position. An undirected g is returned unchanged.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= EdgeList[UndirectedGraph[Graph[{1->2,2->1,2->3}]]]
Out[1]= {1 <-> 2, 2 <-> 3}

In[2]:= EdgeList[UndirectedGraph[Graph[{3->1, 2->3}]]]
Out[2]= {3 <-> 1, 3 <-> 2}
```

### Options (1)

```mathematica
In[3]:= InputForm[UndirectedGraph[Graph[{1,2,3},{1->2,2->1,3->1},EdgeWeight->{2,3,4}]]]
Out[3]= Graph[{1, 2, 3}, {1 <-> 2, 1 <-> 3}, EdgeWeight -> {5, 4}]
```

## Implementation notes

- `Protected`. A non-graph argument is left unevaluated.
- `u -> v` and `v -> u` merge into one edge whose weight is the sum (`Plus`) of theirs.
- Edges are oriented and ordered by `VertexList` position (upper triangle, row-major).
- An undirected `g` is returned unchanged.
- Performance: 1.1–1.6x faster than Mathematica 15 at `10^5` vertices
  (`benchmarks/93-graph-ops-editing`).

**Attributes:** `Protected`.

## References

**See also:** [Plus](../../arithmetic/Plus/), [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

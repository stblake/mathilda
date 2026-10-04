# VertexCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexCount[g] gives the number of vertices in the graph g.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= VertexCount[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[1]= 4

In[2]:= EdgeCount[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[2]= 4

In[3]:= VertexCount[Graph[{}, {}]]
Out[3]= 0

In[4]:= EdgeCount[Graph[{1,2,3},{}]]
Out[4]= 0

In[5]:= EdgeCount[x]
Out[5]= EdgeCount[x]
```

### Applications (3)

K6 has six vertices

```mathematica
In[6]:= VertexCount[CompleteGraph[6]]
Out[6]= 6
```

```mathematica
In[7]:= VertexCount[CycleGraph[4]]
Out[7]= 4
```

Vertices may be arbitrary symbols

```mathematica
In[8]:= VertexCount[Graph[{a, b, c}, {a <-> b}]]
Out[8]= 3
```

## Implementation notes

**Algorithm.** `builtin_vertex_count` returns the number of vertices of `g` as an
integer — a thin reader over the canonical form `Graph[List verts, List edges]`:
the answer is the argument count of the vertex list. When the argument is not a
valid graph it falls through to `hyp_vertex_count`, which handles a `Hypergraph`
(and otherwise returns `NULL`, leaving the expression unevaluated). The companion
`EdgeCount` is the identical reader over the edge list.

**Data structures.** None — it reads `g->data.function.args[0]->...arg_count`
after `graph_is_valid` confirms the shape and populates the memo. The result is a
fresh `Integer`.

**Complexity / limits.** `O(1)` for a canonical graph (plus the one-time
`O(V + E)` validation on first contact with the node). Returns `NULL` for a
non-graph, non-hypergraph argument.

- `Protected`. Cardinalities read from the canonical form; unevaluated on a
  non-graph (see `VertexList`).

**Attributes:** `Protected`.

## References

**See also:** [EdgeCount](../../graphs/EdgeCount/), [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/counts.c`](https://github.com/stblake/mathilda/blob/main/src/graph/counts.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`VertexCount[g]` is the number of vertices, read directly off the canonical
graph form. Vertices can be any expressions, not just integers, and an isolated
vertex (one with no incident edge) still counts.

The parallel reader for edges is `EdgeCount`.

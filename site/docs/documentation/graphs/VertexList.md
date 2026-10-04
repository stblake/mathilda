# VertexList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexList[g] gives the list of vertices of the graph g.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= VertexList[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[1]= {1, 2, 3, 4}

In[2]:= VertexList[Graph[{c->a, a->b}]]
Out[2]= {c, a, b}

In[3]:= VertexList[Graph[{}, {}]]
Out[3]= {}

In[4]:= VertexList[5]
Out[4]= VertexList[5]
```

### Applications (2)

Vertices derived from the edges

```mathematica
In[5]:= VertexList[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]
Out[5]= {1, 2, 3, 4}
```

Vertices may be any expressions

```mathematica
In[6]:= VertexList[PathGraph[{a, b, c}]]
Out[6]= {a, b, c}
```

## Implementation notes

**Algorithm.** `builtin_vertex_list` is a thin reader: for a valid graph it returns a copy of
the graph's vertex `List` — the first argument of the canonical `Graph[List[verts],
List[edges]]` node — in canonical order, with no computation. For a non-graph argument it defers
to the hypergraph reader `hyp_vertex_list`, which handles a `Hypergraph[...]` argument or returns
unevaluated.

**Data structures.** None beyond the deep copy of the stored vertex list; the canonical graph
already holds its vertices in the order `VertexList` reports.

**Complexity / limits.** `O(V)` for the copy. Always exactly one argument.

- `Protected`. The query/representation heads (`VertexList`, `EdgeList`,
  `VertexCount`, `EdgeCount`, `AdjacencyList`, `VertexDegree`,
  `VertexInDegree`, `VertexOutDegree`, `EdgeWeight`) are all thin readers over
  the canonical form and return unevaluated on a non-graph argument.
- Canonical order is the explicit vertex list, or first-appearance order when
  `Graph[e]` derived the vertices from the edges.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/), [VertexCount](../../graphs/VertexCount/), [EdgeCount](../../graphs/EdgeCount/), [AdjacencyList](../../graphs/AdjacencyList/), [VertexDegree](../../graphs/VertexDegree/), [VertexInDegree](../../graphs/VertexInDegree/), [VertexOutDegree](../../graphs/VertexOutDegree/), [EdgeWeight](../../graphs/EdgeWeight/)

- Source: [`src/graph/vertexlist.c`](https://github.com/stblake/mathilda/blob/main/src/graph/vertexlist.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`VertexList[g]` returns the vertices of `g` in canonical order. When a graph is built from edges
alone, the vertices are derived from the edge endpoints in order of first appearance.

Vertices are arbitrary expressions, not just integers; a path on symbols `a`, `b`, `c` lists
them unchanged.

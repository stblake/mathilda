# IncidenceMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IncidenceMatrix[g] gives the vertex-edge incidence matrix of g (oriented: -1 tail, +1 head for directed edges).`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= IncidenceMatrix[PathGraph[3]]
Out[1]= {{1, 0}, {1, 1}, {0, 1}}

In[2]:= IncidenceMatrix[Graph[{1,2,3},{1->2,2->3}]]
Out[2]= {{-1, 0}, {1, -1}, {0, 1}}

In[3]:= IncidenceMatrix[Graph[{1,2,3},{1->2,2<->3}]]
Out[3]= {{-1, 0}, {1, 1}, {0, 1}}

In[4]:= IncidenceMatrix[Graph[{1,2},{}]]
Out[4]= {{}, {}}
```

### Applications (2)

Oriented: -1 tail, +1 head

```mathematica
In[5]:= IncidenceMatrix[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]
Out[5]= {{-1, 0, 1, 0}, {1, -1, 0, 0}, {0, 1, -1, -1}, {0, 0, 0, 1}}
```

```mathematica
In[6]:= MatrixForm[IncidenceMatrix[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]]
Out[6]= MatrixForm[{{-1, 0, 1, 0}, {1, -1, 0, 0}, {0, 1, -1, -1}, {0, 0, 0, 1}}]
```

## Algorithm

incmat.c - IncidenceMatrix[g]: |V| x |E| incidence matrix.

Column j corresponds to edge j (canonical order), row i to vertex i.

```text
  - UndirectedEdge{a,b}: entries (a,j) and (b,j) are 1.
  - DirectedEdge[a,b]:   (a,j) = -1 (tail), (b,j) = 1 (head)  [oriented].
```

Memory (SPEC section 4): returns a freshly-allocated matrix; frees res.

## Implementation notes

**Algorithm.** `builtin_incidence_matrix` builds the dense `|V| x |E|` integer incidence matrix:
row `i` is vertex `i` in `VertexList` order, column `j` is edge `j` in `EdgeList` order. For
each edge it resolves the two endpoint positions and reads the edge kind. An `UndirectedEdge
{a, b}` puts `+1` in both endpoint rows; a `DirectedEdge[a, b]` is oriented, `-1` at the tail and
`+1` at the head. Self-loops never arise, because the `Graph` constructor rejects them. For a
non-graph argument it defers to the hypergraph path `hyp_incidence_matrix`.

**Data structures.** A flat `int` grid (`calloc(n*m)`) filled column by column, then converted
to a `List` of `List`s of boxed integers — a dense matrix, not a `SparseArray`. Endpoint
positions come from a linear `graph_vertex_index` scan.

**Complexity / limits.** `O(V·E)` in both time and space (dense, plus the `O(E·V)` linear
endpoint lookups). The result is a plain nested list suitable for `MatrixForm` or the linear
algebra heads.

- `Protected`. Rows follow `VertexList`, columns follow `EdgeList`. Undirected
  edges mark both endpoints with `1`; directed edges are oriented (`-1` at the
  tail, `+1` at the head).
- Unevaluated on a non-graph.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/), [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/incmat.c`](https://github.com/stblake/mathilda/blob/main/src/graph/incmat.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

The result is a `|V| x |E|` matrix: rows index vertices in `VertexList` order, columns index
edges in `EdgeList` order. A directed edge is oriented — `-1` in its tail row and `+1` in its
head row; an undirected edge puts `+1` in both endpoint rows.

The matrix is a plain nested list, so it feeds straight into `MatrixForm` or the linear algebra
builtins.

# IncidenceMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IncidenceMatrix[g] gives the vertex-edge incidence matrix of g (oriented: -1 tail, +1 head for directed edges).`**

## Examples (4)

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

## Algorithm

incmat.c - IncidenceMatrix[g]: |V| x |E| incidence matrix.

Column j corresponds to edge j (canonical order), row i to vertex i.

```text
  - UndirectedEdge{a,b}: entries (a,j) and (b,j) are 1.
  - DirectedEdge[a,b]:   (a,j) = -1 (tail), (b,j) = 1 (head)  [oriented].
```

Memory (SPEC section 4): returns a freshly-allocated matrix; frees res.

## Implementation notes

- `Protected`. Rows follow `VertexList`, columns follow `EdgeList`. Undirected
  edges mark both endpoints with `1`; directed edges are oriented (`-1` at the
  tail, `+1` at the head).
- Unevaluated on a non-graph.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/), [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

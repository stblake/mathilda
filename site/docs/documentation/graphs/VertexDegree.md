# VertexDegree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexDegree[g] gives the list of vertex degrees; VertexDegree[g,v] gives the degree of vertex v.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= VertexDegree[Graph[{1,2,3},{1<->2,2<->3}]]
Out[1]= {1, 2, 1}

In[2]:= VertexDegree[Graph[{1,2,3},{1<->2,2<->3}], 2]
Out[2]= 2

In[3]:= VertexDegree[Graph[{1,2,3},{1->2,1->3}]]
Out[3]= {2, 1, 1}

In[4]:= VertexDegree[Graph[{1,2},{}]]
Out[4]= {0, 0}
```

## Algorithm

degree.c - VertexDegree, VertexInDegree, VertexOutDegree.

Each accepts VertexDegree[g] (a list of degrees, one per vertex in canonical order) or VertexDegree[g, v] (the degree of a single vertex).

Conventions (documented in docs/spec/builtins/graphs.md):

```text
  - A DirectedEdge[a,b] adds 1 to out(a) and 1 to in(b); total degree of a
    vertex is in + out, so for a purely directed graph total = in + out.
  - An UndirectedEdge[a,b] is incident to both a and b, adding 1 to each of
    their in-, out-, and total degrees (so in = out = total for a purely
    undirected graph).
```

There are no self-loops, so no endpoint is double-counted within one edge.

Memory (SPEC section 4): returns freshly-allocated integers/lists; the evaluator frees res.

## Implementation notes

- `Protected`. Counts every incident edge regardless of direction; for the
  directed split see `VertexInDegree` / `VertexOutDegree`.
- Unevaluated on a non-graph (see `VertexList`) or when `v` is not a vertex.

**Attributes:** `Protected`.

## References

**See also:** [VertexInDegree](../../graphs/VertexInDegree/), [VertexOutDegree](../../graphs/VertexOutDegree/), [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

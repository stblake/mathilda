# VertexOutDegree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexOutDegree[g] / VertexOutDegree[g,v] gives out-degrees (outgoing directed edges; undirected edges count for both).`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= VertexInDegree[Graph[{1,2,3},{1->2,1->3}]]
Out[1]= {0, 1, 1}

In[2]:= VertexOutDegree[Graph[{1,2,3},{1->2,1->3}]]
Out[2]= {2, 0, 0}

In[3]:= VertexInDegree[Graph[{1,2,3},{1->2,1->3}], 2]
Out[3]= 1

In[4]:= VertexOutDegree[Graph[{1,2,3},{1->2,1->3}], 1]
Out[4]= 2

In[5]:= VertexInDegree[Graph[{1,2,3},{1<->2,2->3}]]
Out[5]= {1, 1, 1}

In[6]:= VertexOutDegree[Graph[{1,2,3},{1<->2,2->3}]]
Out[6]= {1, 2, 0}
```

### Applications (3)

Vertex 1 has two out-edges

```mathematica
In[7]:= VertexOutDegree[Graph[{1, 2, 3}, {1 -> 2, 1 -> 3}]]
Out[7]= {2, 0, 0}
```

The out-degree of a single vertex

```mathematica
In[8]:= VertexOutDegree[Graph[{1, 2, 3}, {1 -> 2, 1 -> 3}], 1]
Out[8]= 2
```

An undirected edge counts for out-degree too

```mathematica
In[9]:= VertexOutDegree[CycleGraph[4]]
Out[9]= {2, 2, 2, 2}
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

**Algorithm.** `builtin_vertex_out_degree` gives the out-degree of each vertex
(`VertexOutDegree[g]`, in `VertexList` order) or of a single vertex
(`VertexOutDegree[g, v]`). A `DirectedEdge[a, b]` adds `1` to `out(a)`; an
`UndirectedEdge[a, b]` is incident to both ends and so adds `1` to each of their
in-, out-, and total degrees (hence `in = out = total` for a purely undirected
graph). It shares `degree_dispatch` with `VertexDegree` and `VertexInDegree`,
selecting the `DEG_OUT` accumulator. The whole-graph form makes a single pass over
the edges, pushing each edge's contribution to its endpoints, rather than one
`O(E)` scan per vertex.

**Data structures.** For the list form, a per-vertex `int64 deg[]` accumulator and
a `GraphVIdx` hash index from vertex expression to position, so endpoints resolve
in `O(1)` and a repeated vertex reports the degree of its first occurrence. The
single-vertex form scans the edge list once (`degree_of`).

**Complexity / limits.** `O(V + E)` for the list, `O(E)` for a single vertex —
the earlier per-vertex design was `O(V*E)` (~5 s for 20000 vertices). A
`VertexOutDegree[g, v]` with `v` not a vertex leaves the call unevaluated.

- `Protected`. A `DirectedEdge` adds to the source's out-degree and the target's
  in-degree; an `UndirectedEdge` adds to both the in- and out-degree of each
  endpoint.
- Unevaluated on a non-graph (see `VertexList`).

**Attributes:** `Protected`.

## References

**See also:** [VertexInDegree](../../graphs/VertexInDegree/), [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/degree.c`](https://github.com/stblake/mathilda/blob/main/src/graph/degree.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

`VertexOutDegree[g]` gives the number of edges leaving each vertex, in
`VertexList` order; `VertexOutDegree[g, v]` the value for one vertex. A directed
edge contributes to the out-degree of its tail only; an undirected edge is
incident to both ends, so it contributes to the out-degree of each — which is
why on an undirected graph the out-degree equals the ordinary degree.

The companions are `VertexInDegree` and `VertexDegree` (the total).

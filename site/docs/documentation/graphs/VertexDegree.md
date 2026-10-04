# VertexDegree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexDegree[g] gives the list of vertex degrees; VertexDegree[g,v] gives the degree of vertex v.`**

## Examples (10)

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

### Applications (6)

One degree per vertex, in vertex order

```mathematica
In[5]:= VertexDegree[StarGraph[5]]
Out[5]= {4, 1, 1, 1, 1}
```

The degree of a single vertex

```mathematica
In[6]:= VertexDegree[CycleGraph[6], 1]
Out[6]= 2
```

Every vertex meets the other three

```mathematica
In[7]:= VertexDegree[CompleteGraph[4]]
Out[7]= {3, 3, 3, 3}
```

Directed edges count at both ends

```mathematica
In[8]:= VertexDegree[Graph[{1 -> 2, 1 -> 3, 3 -> 1}]]
Out[8]= {3, 1, 2}
```

The handshake lemma

```mathematica
In[9]:= Total[VertexDegree[PetersenGraph[]]] == 2 EdgeCount[PetersenGraph[]]
Out[9]= True
```

An isolated vertex has degree zero

```mathematica
In[10]:= VertexDegree[Graph[{a, b, c}, {a <-> b}]]
Out[10]= {1, 1, 0}
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

**Algorithm.** `degree_dispatch` serves `VertexDegree[g]` and `VertexDegree[g, v]`, and shares its code with `VertexInDegree` and `VertexOutDegree`. A `DirectedEdge[a, b]` adds 1 to the out-degree of `a` and the in-degree of `b`, and total degree is in plus out. An `UndirectedEdge[a, b]` adds 1 to each endpoint's in, out and total degree. There are no self-loops, so no edge is counted twice at one vertex. A hypergraph argument is handed to `hyp_vertex_degree`.

**Data structures.** The list form makes a single pass over the edges. It resolves endpoints through a `GraphVIdx` hash index into an `int64_t` counter array, then emits a `List` of integers in canonical `VertexList` order. This replaced one `O(E)` scan per vertex, which took about 5 s at 20000 vertices and 40000 edges. The single-vertex form does one `expr_eq` scan over the edge list.

**Complexity / limits.** `O(V + E)` for the full list and `O(E)` for one vertex. The result is not cached. If `v` is not a vertex of `g`, or the argument is not a graph or hypergraph, the call is left unevaluated.

- `Protected`. Counts every incident edge regardless of direction; for the
  directed split see `VertexInDegree` / `VertexOutDegree`.
- Unevaluated on a non-graph (see `VertexList`) or when `v` is not a vertex.

**Attributes:** `Protected`.

## References

**See also:** [VertexInDegree](../../graphs/VertexInDegree/), [VertexOutDegree](../../graphs/VertexOutDegree/), [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/degree.c`](https://github.com/stblake/mathilda/blob/main/src/graph/degree.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`VertexDegree[g]` lists degrees in `VertexList` order, and `VertexDegree[g, v]` gives one vertex. In a directed graph an edge counts at both endpoints, so total degree is in-degree plus out-degree. `VertexInDegree` and `VertexOutDegree` share the same code. A `v` that is not a vertex leaves the call unevaluated.

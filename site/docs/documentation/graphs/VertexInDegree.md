# VertexInDegree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexInDegree[g] / VertexInDegree[g,v] gives in-degrees (incoming directed edges; undirected edges count for both).`**

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

In-degrees in vertex order

```mathematica
In[7]:= VertexInDegree[Graph[{1 -> 2, 1 -> 3, 2 -> 3}]]
Out[7]= {0, 1, 2}
```

The in-degree of a single vertex

```mathematica
In[8]:= VertexInDegree[Graph[{1 -> 2, 1 -> 3, 2 -> 3}], 3]
Out[8]= 2
```

Undirected: in-degree equals total degree

```mathematica
In[9]:= VertexInDegree[CycleGraph[4]]
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

**Algorithm.** `builtin_vertex_in_degree` is `degree_dispatch` in the `DEG_IN`
mode. A `DirectedEdge[a, b]` (i.e. `a -> b`) adds one to the in-degree of `b`
only; an `UndirectedEdge[a, b]` is incident to both ends and contributes one to
the in-degree of each (so in-degree equals out-degree equals total degree for a
purely undirected graph). `VertexInDegree[g, v]` returns the single integer
in-degree of `v` by the per-vertex scan `degree_of`, declining (`NULL`) when `v`
is not a vertex; `VertexInDegree[g]` returns the list of in-degrees in canonical
vertex order. A `Hypergraph` has only a total degree, so the in/out modes return
`NULL` for one.

**Data structures.** The all-vertices form makes a single `O(E)` pass that
*pushes* each edge's contribution onto its endpoints, indexing vertices through a
`GraphVIdx` open-addressing hash — this replaced an earlier pull-style
`O(V·E)` scan (~5 s on a 20000/40000 graph). Counts accumulate in an `int64`
array and are read back through the same index so a repeated vertex reports the
degree of its first occurrence; results are boxed into a `List` of integers. The
single-vertex form uses no auxiliary structure beyond the edge scan.

**Complexity / limits.** `O(V + E)` for `VertexInDegree[g]`, `O(E)` for the
single-vertex query. Self-loops are assumed absent (the constructor forbids
them), so no endpoint is double-counted within one edge. Accepts one or two
arguments; any other arity returns `NULL`.

- `Protected`. A `DirectedEdge` adds to the source's out-degree and the target's
  in-degree; an `UndirectedEdge` adds to both the in- and out-degree of each
  endpoint.
- Unevaluated on a non-graph (see `VertexList`).

**Attributes:** `Protected`.

## References

**See also:** [VertexOutDegree](../../graphs/VertexOutDegree/), [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/degree.c`](https://github.com/stblake/mathilda/blob/main/src/graph/degree.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

A directed edge `a -> b` contributes to the in-degree of `b` only. An undirected
edge is incident to both endpoints, so for a purely undirected graph in-degree
equals out-degree equals total degree.

`VertexInDegree[g]` gives the list of in-degrees in canonical vertex order;
`VertexInDegree[g, v]` gives the single in-degree of `v`, and leaves itself
unevaluated when `v` is not a vertex. A `Hypergraph` has only a total degree, so
the in-degree form does not apply to one.

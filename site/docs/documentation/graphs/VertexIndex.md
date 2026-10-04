# VertexIndex

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexIndex[g, v] gives the position of the vertex v in VertexList[g]; VertexIndex[g, {v1, ...}] gives a list of positions.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= VertexIndex[Graph[{a<->b, b<->c}], b]
Out[1]= 2

In[2]:= VertexIndex[Graph[{a<->b, b<->c}], {c, a}]
Out[2]= {3, 1}

In[3]:= EdgeIndex[CycleGraph[4], 3<->2]
Out[3]= 2

In[4]:= EdgeIndex[CycleGraph[4], 3->4]
Out[4]= 3

In[5]:= EdgeIndex[CycleGraph[4], {1<->2, 4<->1}]
Out[5]= {1, 4}

In[6]:= EdgeIndex[Graph[{1->2}], 2->1]
Out[6]= EdgeIndex[Graph[<2 vertices, 1 edge>], 2 -> 1]
```

### Applications (4)

Position of a vertex in VertexList

```mathematica
In[7]:= VertexIndex[CycleGraph[5], 3]
Out[7]= 3
```

Symbolic vertices work too

```mathematica
In[8]:= VertexIndex[Graph[{a, b, c}, {a <-> b, b <-> c}], c]
Out[8]= 3
```

A list gives a list of positions

```mathematica
In[9]:= VertexIndex[Graph[{a, b, c}, {a <-> b, b <-> c}], {c, a}]
Out[9]= {3, 1}
```

The last vertex of a pentagon

```mathematica
In[10]:= VertexIndex[CycleGraph[5], 5]
Out[10]= 5
```

## Implementation notes

**Algorithm.** `builtin_vertex_index` answers `VertexIndex[g, v]` with the 1-based position of `v` in `VertexList[g]`, and `VertexIndex[g, {v1, v2, ...}]` with the list of positions. It calls `graph_vertex_position`, which probes the memoized vertex hash and returns `-1` for an absent vertex. A single vertex is tried first, so a vertex that is itself a list is still found as a vertex. If any item of a list argument is absent, the whole call is left unevaluated rather than returning a partial list.

**Data structures.** The graph is a `Graph[List, List]` expression tree; the vertex index is a hash from vertex expression to position, built once per graph by `graph_util.c` and reused. The result is plain integers.

**Complexity / limits.** `O(1)` expected per vertex once the index exists (`O(V)` on the first query of a graph). Returns unevaluated for a non-graph, a missing vertex or a wrong argument count.

- `Protected`. A non-graph first argument, or a vertex/edge not in the graph, is
  left unevaluated.
- `EdgeIndex` matches `u -> v` against an undirected edge of an undirected
  graph, as Mathematica does; an undirected edge matches either orientation.

**Attributes:** `Protected`.

## References

**See also:** [EdgeIndex](../../graphs/EdgeIndex/)

- Source: [`src/graph/gops_edit.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_edit.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

The index is 1-based and follows `VertexList`, which need not be sorted. Every vertex in a list argument must exist, otherwise the whole call is left unevaluated.

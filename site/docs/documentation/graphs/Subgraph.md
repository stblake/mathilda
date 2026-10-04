# Subgraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Subgraph[g, {v1, v2, ...}] gives the subgraph of g induced by the listed vertices (elements that are not vertices of g are ignored); Subgraph[g, patt] uses the vertices matching patt. Vertices come in the order given; edges in lower-triangular order of that vertex order. Edge weights are kept.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= EdgeList[Subgraph[Graph[{1,2,3,4},{3<->4,1<->2,2<->3}], {3,2,4}]]
Out[1]= {2 <-> 3, 3 <-> 4}

In[2]:= VertexList[Subgraph[CycleGraph[5], {3,9,1,3,2}]]
Out[2]= {3, 1, 2}

In[3]:= EdgeList[Subgraph[CycleGraph[6], {1,2,3}]]
Out[3]= {1 <-> 2, 2 <-> 3}

In[4]:= EdgeList[Subgraph[CycleGraph[6], _?OddQ]]
Out[4]= {}

In[5]:= Subgraph[CycleGraph[4], {1<->2}]
Out[5]= Subgraph[Graph[<4 vertices, 4 edges>], {TwoWayRule[1, 2]}]
```

### Applications (3)

The induced triangle on three vertices

```mathematica
In[6]:= EdgeList[Subgraph[CompleteGraph[4], {1, 2, 3}]]
Out[6]= {1 <-> 2, 1 <-> 3, 2 <-> 3}
```

Keeps exactly the chosen vertices

```mathematica
In[7]:= VertexCount[Subgraph[CycleGraph[5], {1, 2, 3}]]
Out[7]= 3
```

Only edges with both endpoints kept

```mathematica
In[8]:= EdgeList[Subgraph[CycleGraph[5], {1, 2, 3}]]
Out[8]= {1 <-> 2, 2 <-> 3}
```

## Implementation notes

**Algorithm.** `builtin_subgraph` gives the subgraph of `g` **induced** by a set
of vertices: every edge of `g` with both endpoints in the set is kept. The second
argument may be an explicit vertex list (elements that are not vertices of `g`
are ignored; an edge element makes it decline, since edge-induced subgraphs are
not handled), a pattern (every matching vertex is selected via `gops_matchq`), or
a single vertex. Selected vertices are given result positions in the order
listed, then `build_induced` keeps each edge whose endpoints both survive,
emitting edges in lower-triangular adjacency order of that vertex order. Edge
weights are kept.

**Data structures.** A `pos[]` map from original to result vertex index, a
`GopsView` of the graph, and an induced-subgraph builder that walks each kept
vertex's incidence, carrying the aligned `EdgeWeight` entries and sharing the
edge nodes. The result is seeded into the memo (`gops_graph_new`).

**Complexity / limits.** `O(V + E)` plus one position lookup per selected vertex.
A single-vertex argument that is not a vertex of `g`, or an edge-list argument,
leaves the call unevaluated. The companion `NeighborhoodGraph` induces on a
distance ball instead of an explicit set.

- `Protected`. A non-graph first argument is left unevaluated.
- Vertices appear in the given order; non-vertices are ignored and repeats dropped.
- Edges are emitted at their later endpoint in that order, following each
  vertex's incidence order (out-edges and undirected edges, then in-edges).
- Weights are kept.
- An edge list (the edge-induced subgraph form) is not supported; the call is
  left unevaluated.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_edit.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_edit.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`Subgraph[g, vs]` is the *vertex-induced* subgraph: it keeps the listed vertices
and exactly those edges of `g` whose endpoints both lie in the list. In a cycle,
inducing on three consecutive vertices keeps the two edges between them but drops
the wrap-around edge.

Vertices come out in the order given; list elements that are not vertices of `g`
are ignored. Edge weights are preserved.

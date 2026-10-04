# NeighborhoodGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NeighborhoodGraph[g, v] gives the subgraph of g induced by v and its neighbours; NeighborhoodGraph[g, v, k] by the vertices within distance k of v (k a non-negative integer or Infinity), edge direction ignored. NeighborhoodGraph[g, {v1, ...}, k] uses several centres. Vertices: the centres, then each centre's new vertices in VertexList order.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= VertexList[NeighborhoodGraph[PathGraph[Range[6]], {1,6}]]
Out[1]= {1, 6, 2, 5}

In[2]:= VertexList[NeighborhoodGraph[PathGraph[Range[6]], 3, 2]]
Out[2]= {3, 1, 2, 4, 5}

In[3]:= EdgeList[NeighborhoodGraph[Graph[{1->2,3->1,2->4}], 1]]
Out[3]= {1 -> 2, 3 -> 1}

In[4]:= VertexList[NeighborhoodGraph[PathGraph[Range[6]], 2, Infinity]]
Out[4]= {2, 1, 3, 4, 5, 6}

In[5]:= VertexList[NeighborhoodGraph[PathGraph[Range[6]], 9]]
Out[5]= {}
```

### Applications (4)

Vertex 1 and its two cycle neighbours

```mathematica
In[6]:= VertexList[NeighborhoodGraph[CycleGraph[6], 1]]
Out[6]= {1, 2, 6}
```

The ball of radius two around vertex 1

```mathematica
In[7]:= VertexList[NeighborhoodGraph[CycleGraph[6], 1, 2]]
Out[7]= {1, 2, 3, 5, 6}
```

The centre first, then its ball in order

```mathematica
In[8]:= VertexList[NeighborhoodGraph[PathGraph[{1, 2, 3, 4, 5}], 3, 1]]
Out[8]= {3, 2, 4}
```

In K5 one step reaches every vertex

```mathematica
In[9]:= EdgeCount[NeighborhoodGraph[CompleteGraph[5], 1]]
Out[9]= 10
```

## Implementation notes

**Algorithm.** `builtin_neighborhood_graph` returns the subgraph of `g` induced by
a set of centres together with every vertex within graph distance `k` of them
(default `k = 1`; `k = Infinity` and any non-negative integer are accepted, edge
**direction ignored**). The centre argument is a single vertex or a list of them;
non-vertices are silently dropped. From each centre it runs a BFS over the
all-direction incidence (`gops_inc_build` with `GOPS_INC_ALL`), stamped per centre
so the balls of several centres do not interfere, stopping when a vertex's distance
reaches `k`. Result vertices are listed centres-first, then each centre's newly
reached ball in `VertexList` order — a `qsort` of the ball, or a linear rescan over
all vertices when the ball exceeds roughly `1/8` of them. `build_induced` then keeps
every edge of `g` whose both endpoints were selected.

**Data structures.** A `GopsView` exposes the memo's endpoint arrays (`eu`/`ev`/
`edir`); `gops_inc_build` turns them into a CSR incidence (`GopsInc`, neighbours of
each vertex in one `start`/`nbr` pair). The BFS uses per-centre `mark` (a stamp),
`dist`, `queue`, and `ball` arrays plus a `pos[]` map from graph vertex to result
position. `build_induced` emits each kept edge once — at its later endpoint, in
lower-triangular adjacency order — reusing the shared edge, weight, and vertex nodes
through `expr_copy`, and carries an `EdgeWeight` list aligned to the surviving edges.

**Complexity / limits.** `O(V + E)` per centre for the BFS and `O(V + E)` to induce
the subgraph; each edge is emitted once. `k = Infinity` yields the whole connected
hull of the centres. Arity other than two or three, or a `k` that is neither a
non-negative integer nor `Infinity`, returns `NULL`. `res` is borrowed and never
modified.

- `Protected`. A non-graph first argument is left unevaluated.
- Distance ignores edge direction.
- Vertices: the centres, then each centre's new vertices in `VertexList` order;
  edges are ordered as for `Subgraph`.
- Non-vertex centres are ignored.
- `k = Infinity` is accepted (Mathematica leaves it unevaluated).

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/), [Subgraph](../../graphs/Subgraph/)

- Source: [`src/graph/gops_edit.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_edit.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`NeighborhoodGraph[g, v, k]` is the subgraph of `g` induced by `v` and every vertex
within graph distance `k` of it (`k` defaults to `1`, and `Infinity` is allowed).
Edge direction is ignored when measuring distance. The centre may instead be a list
of vertices, and non-vertices among them are quietly skipped. The result is an
opaque `Graph`, read through `VertexList`, `EdgeList`, `EdgeCount`, and the others.

Vertices come out centres-first, then each centre's newly reached ball in
`VertexList` order — so `NeighborhoodGraph[CycleGraph[6], 1]` lists `1` before its
neighbours `2` and `6`. With `k` large enough the result is the whole connected
component of the centres.

# AdjacencyList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AdjacencyList[g] gives the adjacency list of g; AdjacencyList[g,v] gives the vertices adjacent to v (successors for directed edges).`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= AdjacencyList[Graph[{1,2,3},{1<->2,2<->3}]]
Out[1]= {{2}, {1, 3}, {2}}

In[2]:= AdjacencyList[Graph[{1,2,3},{1<->2,2<->3}], 2]
Out[2]= {1, 3}

In[3]:= AdjacencyList[Graph[{1,2,3},{1->2,3->1}]]
Out[3]= {{2}, {}, {1}}

In[4]:= AdjacencyList[Graph[{1,2,3},{1->2,3->1}], 1]
Out[4]= {2}

In[5]:= AdjacencyList[5]
Out[5]= AdjacencyList[5]
```

### Applications (4)

Each vertex of a 4-cycle has two neighbours

```mathematica
In[6]:= AdjacencyList[CycleGraph[4]]
Out[6]= {{2, 4}, {1, 3}, {2, 4}, {3, 1}}
```

The neighbours of a single vertex

```mathematica
In[7]:= AdjacencyList[CompleteGraph[4], 1]
Out[7]= {2, 3, 4}
```

Directed: only successors count

```mathematica
In[8]:= AdjacencyList[Graph[{1 -> 2, 1 -> 3, 2 -> 3}], 1]
Out[8]= {2, 3}
```

A directed sink has no successors

```mathematica
In[9]:= AdjacencyList[Graph[{1 -> 2, 2 -> 3}], 3]
Out[9]= {}
```

## Algorithm

adjlist.c - AdjacencyList[g] and AdjacencyList[g, v].

Neighbors of v: successors for directed edges (v -> u yields u), and both endpoints for undirected edges. This matches the row convention of AdjacencyMatrix (a 1 in row v, column u means v is adjacent to u). Neighbors are returned in first-appearance order, de-duplicated.

```text
AdjacencyList[g]     -> {neighbors(v1), neighbors(v2), ...} in vertex order.
AdjacencyList[g, v]  -> neighbors(v).
```

Memory (SPEC section 4): returns freshly-allocated lists; evaluator frees res.

## Implementation notes

**Algorithm.** `builtin_adjacency_list` returns the neighbours of each vertex, following
the same orientation convention as `AdjacencyMatrix`: for a `DirectedEdge[v, u]` only the
successor `u` is a neighbour of `v`, while an `UndirectedEdge` contributes both endpoints.
`AdjacencyList[g]` gives `{neighbours(v1), neighbours(v2), ...}` in canonical vertex order;
`AdjacencyList[g, v]` gives just the neighbours of `v`. For each vertex the helper
`neighbors_of` scans `g`'s edge list once, comparing endpoints by `expr_eq` and appending
each new neighbour in **first-appearance order**, de-duplicated (so a vertex reachable by
both an undirected and a redundant directed edge is listed once).

**Data structures.** The graph is the canonical `Graph[List verts, List edges]` `Expr`
tree, read directly — `neighbors_of` allocates a scratch `Expr*` buffer of size `2·|E|`,
`expr_copy`s each surviving neighbour into a fresh `List`, and the one-argument form wraps
the per-vertex lists in an outer `List`. Vertex membership for the two-argument form is
checked with `graph_vertex_index` (a linear `expr_eq` scan over `verts`).

**Complexity / limits.** `neighbors_of` is `O(E)` per vertex with an inner dedup scan, so
`AdjacencyList[g]` is `O(V·E)` in the worst case (it reads the raw edge list rather than
the shared `GraphAdj` CSR). The head returns `NULL` (unevaluated) when `g` is not a valid
graph, and `AdjacencyList[g, v]` returns `NULL` when `v` is not a vertex of `g`.

- `Protected`. Directed edges contribute successors (`v -> u` makes `u` a
  neighbor of `v`, not the reverse); undirected edges go both ways.
- Unevaluated on a non-graph (see `VertexList`) or when `v` is not a vertex.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/adjlist.c`](https://github.com/stblake/mathilda/blob/main/src/graph/adjlist.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

Neighbours follow edge orientation: for a directed edge `v -> u` the successor `u` is a
neighbour of `v` (but `v` is not a neighbour of `u`), while an undirected edge makes each
endpoint a neighbour of the other. This matches the row convention of `AdjacencyMatrix`.

Neighbours are returned in first-appearance order and de-duplicated. `AdjacencyList[g]`
lists all vertices' neighbourhoods in `VertexList` order; `AdjacencyList[g, v]` gives just
one. A non-graph argument, or a `v` that is not a vertex of `g`, leaves the expression
unevaluated.

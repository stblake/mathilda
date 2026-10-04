# PetersenGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PetersenGraph[] gives the Petersen graph; PetersenGraph[n, k] the generalized Petersen graph: inner vertices 1..n with i joined to i+k (mod n), outer cycle n+1..2n, and spokes i to n+i.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= PetersenGraph[]
Out[1]= Graph[<10 vertices, 15 edges>]

In[2]:= VertexDegree[PetersenGraph[]]
Out[2]= {3, 3, 3, 3, 3, 3, 3, 3, 3, 3}

In[3]:= GraphDiameter[PetersenGraph[]]
Out[3]= 2

In[4]:= EdgeList[PetersenGraph[4, 1]]
Out[4]= {1 <-> 2, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 6, 3 <-> 4, 3 <-> 7, 4 <-> 8, 5 <-> 6, 5 <-> 8, 6 <-> 7, 7 <-> 8}
```

### Applications (5)

The Petersen graph has 10 vertices

```mathematica
In[5]:= VertexCount[PetersenGraph[]]
Out[5]= 10
```

And 15 edges

```mathematica
In[6]:= EdgeCount[PetersenGraph[]]
Out[6]= 15
```

It is cubic: every vertex has degree 3

```mathematica
In[7]:= VertexDegree[PetersenGraph[]]
Out[7]= {3, 3, 3, 3, 3, 3, 3, 3, 3, 3}
```

The generalized form GP(n, k) has 2n vertices labelled 1..2n

```mathematica
In[8]:= VertexList[PetersenGraph[3, 1]]
Out[8]= {1, 2, 3, 4, 5, 6}
```

A larger generalized Petersen graph

```mathematica
In[9]:= VertexCount[PetersenGraph[7, 2]]
Out[9]= 14
```

## Implementation notes

**Algorithm.** `builtin_petersen_graph` emits the generalized Petersen graph
`GP(n, k)` on `2n` vertices. `PetersenGraph[]` is the default `GP(5, 2)` — the
Petersen graph itself. For `PetersenGraph[n, k]` the three edge families are added
with `pairs_add` over `i = 0 … n−1`: the inner "star" `i — (i+k) mod n`, the outer
cycle `(n+i) — (n + (i+1) mod n)`, and the spokes `i — (n+i)`. The result is handed
to `pairs_graph`, which renumbers the vertices to `1 … 2n` (inner `1 … n`, outer
`n+1 … 2n`) and deduplicates, so the doubled inner/outer edges at `n = 2` collapse
to a simple graph.

**Data structures.** A `Pairs` endpoint-pair accumulator filled in one pass, then
converted to the `Graph` expression by `pairs_graph` (which sorts and dedups the
edge list).

**Complexity / limits.** `O(n)` edges emitted and a linear build. The arguments
must satisfy `n ≥ 0`, `k ≥ 0`, `k` not a multiple of `n` (`k mod n ≠ 0`, so the
inner edges form a genuine circulant and not self-loops), and `n ≤
GEN_MAX_VERTICES/2`; otherwise the call is left unevaluated. The standard Petersen
graph `GP(5, 2)` is 3-regular on 10 vertices with 15 edges.

- Undirected on `1..2n`: the inner star is `1..n`, the outer cycle `n+1..2n`;
  the edge list is the sorted list of pairs `{i, j}`, `i < j`, identical to
  Wolfram's `EdgeList`.
- `PetersenGraph[]` is `PetersenGraph[5, 2]`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- M. E. Watkins, *A theorem on Tait colorings with an application to the generalized Petersen graphs*, J. Combin. Theory **6** (1969) 152-164.
- Source: [`src/graph/gmet_generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graphplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphplot.c)

## Notes & additional examples

### Notes

`PetersenGraph[]` is the classic Petersen graph `GP(5, 2)`: an inner pentagram
joined by spokes to an outer pentagon, 3-regular on 10 vertices. `PetersenGraph[n,
k]` gives the generalized Petersen graph `GP(n, k)` on `2n` vertices, with inner
vertices `1 … n` joined at step `k` and outer vertices `n+1 … 2n` forming a cycle.

The construction requires that `k` not be a multiple of `n`, so the inner edges are
real chords rather than self-loops. The result is an opaque `Graph` object; query a
property (`VertexCount`, `EdgeCount`, `VertexDegree`, `VertexList`) to inspect it.

# GraphTriangleCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphTriangleCount[g] gives the number of triangles of an undirected graph g, or of directed 3-cycles of a directed graph. O(m^1.5) degree-ordered triangle listing.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= GraphTriangleCount[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[1]= 1

In[2]:= GraphTriangleCount[CompleteGraph[5]]
Out[2]= 10

In[3]:= GraphTriangleCount[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= 1

In[4]:= GraphTriangleCount[Graph[{1->2,2->3,1->3}]]
Out[4]= 0

In[5]:= GraphTriangleCount[Graph[{1<->2, 2->3}]]
Out[5]= GraphTriangleCount[Graph[<3 vertices, 2 edges>]]
```

### Applications (4)

C(4,3) = 4 triangles

```mathematica
In[6]:= GraphTriangleCount[CompleteGraph[4]]
Out[6]= 4
```

C(5,3) = 10 triangles

```mathematica
In[7]:= GraphTriangleCount[CompleteGraph[5]]
Out[7]= 10
```

A cycle has none

```mathematica
In[8]:= GraphTriangleCount[CycleGraph[5]]
Out[8]= 0
```

The Petersen graph is triangle-free

```mathematica
In[9]:= GraphTriangleCount[PetersenGraph[]]
Out[9]= 0
```

## Implementation notes

**Algorithm.** `builtin_graph_triangle_count` gives the number of triangles of an
undirected graph, or of directed 3-cycles `u->v->w->u` of a directed graph. The
core `tri_compute` lists triangles over the underlying simple graph with the
**degree-ordered orientation**: each edge is oriented from the endpoint of lower
`(degree, index)` rank to the higher, which bounds every oriented out-degree by
`O(sqrt m)` and the total listing work by `O(m^1.5)` (Chiba-Nishizeki / Latapy).
For each oriented edge `u -> v` it intersects `u`'s and `v`'s forward
neighbourhoods via a mark array. For a directed graph each oriented edge carries
a two-bit record of which of the two arcs exist, so the 3-cycle test is two
bit-ANDs per candidate; a transitive triple counts `0` and a doubly-linked
triangle counts `2`.

**Data structures.** CSR oriented adjacency built from the memo's integer
endpoints: `ooff`/`oadj`, plus a per-edge flag byte `ofl` for the directed case.
Rows are processed on the thread team with per-thread mark arrays and counters,
and the undirected build uses index-order orientation when `maxdeg^2 <= 4m`
(where the two bounds agree) to avoid the degree sort.

**Complexity / limits.** `O(m^1.5)` time, exact `Integer` result. `EdgeWeight` is
ignored; a mixed graph is left unevaluated.

- Exact, weights ignored, mixed graphs unevaluated (as in Wolfram).
- Undirected: the number of 3-cliques. Directed (*reverse-engineered*):
  triangles are directed 3-cycles, so a transitive triple `1->2, 2->3, 1->3`
  does not count.
- **Algorithm:** triangle listing with an acyclic orientation (degree order,
  or index order when `maxdeg² <= 4m`), O(m^1.5) worst case; directed 3-cycles
  via two-bit arc flags per edge.

**Attributes:** `Protected`.

## References

- N. Chiba and T. Nishizeki, *Arboricity and subgraph listing algorithms*, SIAM J. Comput. **14** (1985) 210-223.
- M. Latapy, *Main-memory triangle computations for very large (sparse (power-law)) graphs*, Theoret. Comput. Sci. **407** (2008) 458-473.
- Source: [`src/graph/gmet_cluster.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_cluster.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

`GraphTriangleCount[g]` counts triangles (undirected) or directed 3-cycles
`u->v->w->u` (directed). A complete graph `K_n` has `C(n, 3)` triangles.

The count is exact and uses a degree-ordered triangle listing that runs in
`O(m^1.5)`. `EdgeWeight` is ignored; a mixed graph is left unevaluated.

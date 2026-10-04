# GridGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GridGraph[{m, n}] gives the m x n grid graph; GridGraph[{n1, ..., nk}] the k-dimensional grid. The first coordinate varies fastest in the vertex numbering.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= GridGraph[{3, 2}]
Out[1]= Graph[<6 vertices, 7 edges>]

In[2]:= EdgeList[GridGraph[{3, 2}]]
Out[2]= {1 <-> 2, 1 <-> 4, 2 <-> 3, 2 <-> 5, 3 <-> 6, 4 <-> 5, 5 <-> 6}

In[3]:= GridGraph[{2, 2, 2}]
Out[3]= Graph[<8 vertices, 12 edges>]

In[4]:= GridGraph[{5}]
Out[4]= Graph[<5 vertices, 4 edges>]

In[5]:= GridGraph[Table[2, {26}]]
Out[5]= GridGraph[{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}]
```

### Applications (6)

A 2 by 3 grid, first coordinate varying fastest

```mathematica
In[6]:= EdgeList[GridGraph[{2, 3}]]
Out[6]= {1 <-> 2, 1 <-> 3, 2 <-> 4, 3 <-> 4, 3 <-> 5, 4 <-> 6, 5 <-> 6}
```

A bare integer gives a path on 4 vertices

```mathematica
In[7]:= EdgeList[GridGraph[4]]
Out[7]= {1 <-> 2, 2 <-> 3, 3 <-> 4}
```

Corners 2, edge midpoints 3, centre 4

```mathematica
In[8]:= VertexDegree[GridGraph[{3, 3}]]
Out[8]= {2, 3, 2, 3, 4, 3, 2, 3, 2}
```

A three-dimensional grid has the product of the dimensions

```mathematica
In[9]:= VertexCount[GridGraph[{2, 3, 4}]]
Out[9]= 24
```

The sum over axes of (n_i - 1) times the product of the others

```mathematica
In[10]:= EdgeCount[GridGraph[{3, 4, 5}]]
Out[10]= 133
```

The cube graph, identical to HypercubeGraph[3]

```mathematica
In[11]:= EdgeList[GridGraph[{2, 2, 2}]]
Out[11]= {1 <-> 2, 1 <-> 3, 1 <-> 5, 2 <-> 4, 2 <-> 6, 3 <-> 4, 3 <-> 7, 4 <-> 8, 5 <-> 6, 5 <-> 7, 6 <-> 8, 7 <-> 8}
```

## Implementation notes

**Algorithm.** `builtin_grid_graph` accepts one argument: an integer `n` (a path of `n` vertices) or a list `{n1, ..., nk}` of up to 32 positive dimensions. Vertices are numbered `1 + x1 + n1*x2 + n1*n2*x3 + ...`, so the first coordinate varies fastest. For each vertex and each axis `i` with room to step (`x_i + 1 < n_i`), it emits the edge to the vertex one stride further along that axis.

**Data structures.** Per-axis strides are held in small `int64_t[32]` arrays and edges go into the packed-64-bit `Pairs` buffer, which is already in sorted order for this loop. `pairs_graph` then returns the canonical `Graph[Range[N], {UndirectedEdge[u, v], ...}]` tree.

**Complexity / limits.** `O(N * k)` for `N = n1...nk` vertices; the exact edge count `sum_i (n_i - 1) N / n_i` is computed first so an over-cap grid is refused before any allocation. More than 10^8 vertices, 5 x 10^7 edges or 32 dimensions, and any non-positive or non-integer dimension, leave the call unevaluated.

- Undirected on `1..N`, vertices numbered with the first coordinate varying
  fastest; the edge list is the sorted list of pairs `{i, j}`, `i < j`,
  identical to Wolfram's `EdgeList`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limits (shared by every graph family of this module): a family
  with more than 10^8 vertices or 5×10^7 edges is left unevaluated rather than
  allocating gigabytes (e.g. `GridGraph[Table[2, {26}]]`, ~8.7×10^8 edges).

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graphplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphplot.c)

## Notes & additional examples

### Notes

`GridGraph[{n1, ..., nk}]` is the Cartesian product of paths of lengths `n1, ..., nk`. Vertex `1 + x1 + n1 x2 + n1 n2 x3 + ...` sits at coordinates `(x1, x2, ...)`, so the first coordinate runs fastest.

Dimensions must be positive integers (up to 32 of them). Very large grids, over 10^8 vertices or 5 x 10^7 edges, are refused and left unevaluated.

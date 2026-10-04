# HypergraphDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphDistance[h, u, v] gives the least number of hyperedges in a chain joining vertices u and v, or Infinity. HypergraphDistance[h, u] gives the distances from u to every vertex, in VertexList order.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= HyperedgeDistance[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 1, 3]
Out[1]= 2

In[2]:= HyperedgeDistance[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 1]
Out[2]= {0, 1, 2, Infinity}

In[3]:= HyperedgeDistance[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 1, All, 2]
Out[3]= {0, 1, 2, Infinity}

In[4]:= HypergraphDistance[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 1, 6]
Out[4]= 3

In[5]:= HypergraphDistance[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 1]
Out[5]= {0, 1, 1, 2, 3, 3, Infinity}

In[6]:= HyperedgeDistance[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 1, 9]
Out[6]= HyperedgeDistance[Hypergraph[<7 vertices, 4 hyperedges>], 1, 9]
```

### Applications (4)

```mathematica
In[7]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
Out[7]= Hypergraph[<7 vertices, 4 hyperedges>]
```

Fewest hyperedge hops from vertex 1 to vertex 6

```mathematica
In[8]:= HypergraphDistance[h, 1, 6]
Out[8]= 3
```

Distances from vertex 1 to every vertex, in VertexList order

```mathematica
In[9]:= HypergraphDistance[h, 1]
Out[9]= {0, 1, 1, 2, 3, 3, Infinity}
```

The distance from a vertex to itself is 0

```mathematica
In[10]:= HypergraphDistance[h, 1, 1]
Out[10]= 0
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_distance` gives the least number of hyperedges
in a chain of vertices `u = x_0, x_1, ..., x_k = v` with consecutive vertices
sharing a hyperedge — the shortest-path distance in the clique expansion
(`HypergraphCliqueExpansion`), `Infinity` when none, `0` for `u == v`. It runs a
breadth-first search over vertices that, instead of materialising the clique
graph, expands each hyperedge at most once (an `eexp` stamp): popping a vertex
walks its incidence list, and each not-yet-expanded hyperedge relaxes all its
distinct members. `HypergraphDistance[h, u]` returns distances from `u` to every
vertex in `VertexList` order.

**Data structures.** The distinct-vertex CSR `soff/sv` and incidence CSR
`voff/ve`; a BFS queue, a vertex distance array, and an `eexp` hyperedge-expanded
bitmap; `dist_list`/`dist_expr` render the result (packed `int64` when all finite,
`Infinity` symbols otherwise).

**Complexity / limits.** `O(Σ|e|)` — each hyperedge's members are scanned once.
An unknown source or target vertex leaves the call unevaluated.

- `Infinity` when no walk exists; `0` for `i == j` (resp. `u == v`).
- `HyperedgeDistance` runs a BFS directly on the incidence structure — the line
  graph is never materialised. With `s = 1` it is linear: each vertex's
  incidence list is expanded at most once per BFS; with `s ≥ 2` it is
  `O(Σ_v deg(v)²)`. The all-targets form is packed when all distances are
  finite.
- `HypergraphDistance` is the distance in the clique expansion
  (`HypergraphCliqueExpansion`).
- An out-of-range hyperedge index or an unknown vertex leaves the call
  unevaluated.
- Benchmark (experiment 96), single source on 10⁵ hyperedges:
  `HyperedgeDistance` 2.8 ms warm / 3.3 ms cold (Mathematica 7255 ms, xgi
  1048 ms); `HypergraphDistance` 2.9 ms warm / 3.5 ms cold (Mathematica
  5791 ms, xgi 1057 ms).

**Attributes:** `Protected`.

## References

**See also:** [HyperedgeDistance](../../hypergraphs/HyperedgeDistance/), [HypergraphCliqueExpansion](../../hypergraphs/HypergraphCliqueExpansion/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`HypergraphDistance[h, u, v]` is the least number of hyperedges in a chain of
vertices `u = x_0, x_1, ..., x_k = v` in which consecutive vertices share a
hyperedge — equivalently the shortest-path distance in the clique expansion
`HypergraphCliqueExpansion[h]`. It is `0` for `u == v` and `Infinity` when `u` and
`v` lie in different components. `HypergraphDistance[h, u]` gives the whole
distance vector from `u`.

The BFS walks vertices but expands each hyperedge at most once, so it never
materialises the clique graph and stays linear in the total incidence. Contrast
`HyperedgeDistance`, which measures distance *between hyperedges*. An unknown
vertex leaves the call unevaluated.

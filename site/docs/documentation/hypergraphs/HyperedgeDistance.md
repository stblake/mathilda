# HyperedgeDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HyperedgeDistance[h, i, j] gives the length of a shortest walk of intersecting hyperedges from hyperedge i to hyperedge j, or Infinity. HyperedgeDistance[h, i, j, s] uses s-walks (consecutive hyperedges share at least s vertices). HyperedgeDistance[h, i] or [h, i, All, s] gives the distances from i to every hyperedge.`**

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

Shortest walk of intersecting hyperedges, 1 -> 2 -> 3

```mathematica
In[8]:= HyperedgeDistance[h, 1, 3]
Out[8]= 2
```

Distances from hyperedge 1 to every hyperedge

```mathematica
In[9]:= HyperedgeDistance[h, 1]
Out[9]= {0, 1, 2, Infinity}
```

Hyperedge 4 = {7} is unreachable

```mathematica
In[10]:= HyperedgeDistance[h, 1, 4]
Out[10]= Infinity
```

## Implementation notes

**Algorithm.** `builtin_hyperedge_distance` gives the length of a shortest walk of
intersecting hyperedges from `i` to `j` (an s-walk, consecutive hyperedges sharing
at least `s` vertices, when a fourth argument `s` is given); `Infinity` if none,
`0` for `i == j`. `edge_bfs` runs a breadth-first search directly on the incidence
structure — the line graph is never materialised. For `s = 1` it expands each
vertex's incidence list at most once per BFS (a `vexp` stamp), which keeps it
linear; for `s >= 2` it expands neighbours through `s_neighbours`. The one-source
forms `HyperedgeDistance[h, i]` and `[h, i, All, s]` return the whole distance
vector (stopping early only when a specific target is named).

**Data structures.** The distinct-vertex CSR `soff/sv` and incidence CSR
`voff/ve`; a BFS queue and a `vexp` vertex-expanded bitmap (`s = 1`) or
`cnt`/`touched`/`nb` overlap scratch (`s >= 2`); `dist_list` packs the all-targets
answer into an `int64` buffer when every distance is finite, else emits a `List`
with `Infinity` entries.

**Complexity / limits.** `s = 1` is `O(Σ|e|)`; `s >= 2` is `O(Σ_v deg(v)²)`. An
out-of-range hyperedge index leaves the call unevaluated; `s` must be a positive
integer.

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

**See also:** [HypergraphDistance](../../hypergraphs/HypergraphDistance/), [HypergraphCliqueExpansion](../../hypergraphs/HypergraphCliqueExpansion/)

- S. G. Aksoy, C. Joslyn, C. Ortiz Marrero, B. Praggastis and E. Purvine, *Hypernetwork science via high-order hypergraph walks*, EPJ Data Science **9**:16 (2020).
- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`HyperedgeDistance[h, i, j]` is the number of steps in a shortest walk of
intersecting hyperedges from `i` to `j`, `0` when `i == j` and `Infinity` when no
walk exists. `HyperedgeDistance[h, i, j, s]` uses s-walks (consecutive hyperedges
sharing at least `s` vertices); `HyperedgeDistance[h, i]` and
`HyperedgeDistance[h, i, All, s]` give the distances from `i` to every hyperedge.

The search is a BFS run directly on the incidence structure — the line graph is
never built — so the single-source form is linear in the total incidence for
`s = 1`. Hyperedge indices are 1-based `EdgeList` positions; an out-of-range index
leaves the call unevaluated.

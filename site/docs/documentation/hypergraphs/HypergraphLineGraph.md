# HypergraphLineGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphLineGraph[h] gives the line graph of h: vertices 1..m, i<->j when hyperedges i and j intersect. HypergraphLineGraph[h, s] gives the s-line graph, joining hyperedges that share at least s vertices.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[HypergraphLineGraph[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[1]= {1 <-> 2, 2 <-> 3}

In[2]:= EdgeList[HypergraphLineGraph[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 2]]
Out[2]= {1 <-> 2, 2 <-> 3}

In[3]:= InputForm[HypergraphLineGraph[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 3]]
Out[3]= Graph[{1, 2, 3, 4}, {}]

In[4]:= HypergraphLineGraph[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 0]
Out[4]= HypergraphLineGraph[Hypergraph[<5 vertices, 4 hyperedges>], 0]
```

### Applications (3)

```mathematica
In[5]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
Out[5]= Hypergraph[<7 vertices, 4 hyperedges>]
```

An edge joins hyperedges i and j when they intersect

```mathematica
In[6]:= EdgeList[HypergraphLineGraph[h]]
Out[6]= {1 <-> 2, 2 <-> 3}
```

One vertex per hyperedge, even an isolated one

```mathematica
In[7]:= VertexCount[HypergraphLineGraph[h]]
Out[7]= 4
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_line_graph` gives the line graph on the
hyperedge indices `1..m`, with `i <-> j` when hyperedges `i` and `j` intersect; a
second argument `s` requests the **s-line graph** (`|e_i ∩ e_j| >= s`, after Aksoy
et al.). Intersections are counted through the shared `s_neighbours` helper: for
hyperedge `i`, it walks `i`'s distinct vertices and, through the vertex→hyperedge
incidence CSR, increments a per-hyperedge counter, collecting each touched `j > i`
whose overlap reaches `s`, then zeroes the scratch. Edges are emitted in `(i, j)`
order (neighbours sorted). Every hyperedge index is a vertex, even when isolated.

**Data structures.** The distinct-vertex CSR `soff/sv` and the incidence CSR
`voff/ve` (`hyp_view(..., 1)`); reusable `cnt`/`touched`/`nb` scratch arrays left
zeroed between hyperedges; an `EVec` of `UndirectedEdge` expressions. The result
is a `Graph` on the index integers `1..m`.

**Complexity / limits.** `O(Σ_v deg(v)²)` — the number of hyperedge pairs meeting
at a vertex, which bounds the 1-line graph's edge count anyway. `s` must be a
positive integer, else the call is left unevaluated.

- s-line graphs follow Aksoy, Joslyn, Ortiz Marrero, Praggastis, Purvine,
  "Hypernetwork science via high-order hypergraph walks" (2020). Hyperedges
  are read as sets.
- Edges ordered by `i`, then `j`. Every hyperedge index is a vertex, even when
  isolated.
- `O(Σ_v deg(v)²)` — the number of hyperedge pairs meeting at a vertex, which
  bounds the 1-line graph's size anyway.
- `s` must be a positive integer; otherwise unevaluated.
- Benchmark (experiment 96): 10⁵ hyperedges 47.0 ms warm / 48.5 ms cold
  (Mathematica 819 ms, xgi 1117 ms).

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/)

- S. G. Aksoy, C. Joslyn, C. Ortiz Marrero, B. Praggastis and E. Purvine, *Hypernetwork science via high-order hypergraph walks*, EPJ Data Science **9**:16 (2020).
- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

The line graph puts one vertex per hyperedge (`1..m`) and joins `i` and `j` when
hyperedges `i` and `j` intersect. `HypergraphLineGraph[h, s]` gives the **s-line
graph** of Aksoy et al.: `i <-> j` only when the two hyperedges share at least `s`
vertices, the overlap notion that underlies high-order hypergraph walks.

Hyperedges are read as sets, edges are listed in `(i, j)` order, and every
hyperedge index is a vertex even when isolated. The overlaps are counted straight
through the incidence lists — the full pairwise intersection matrix is never
formed — so the cost is the number of hyperedge pairs meeting at a vertex. `s` must
be a positive integer.

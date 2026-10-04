# HyperedgeConnectedComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HyperedgeConnectedComponents[h] gives the components of the hyperedges of h under intersection, as Lists of 1-based hyperedge indices. HyperedgeConnectedComponents[h, s] gives the s-connected components: hyperedges with at least s vertices, joined when they share at least s.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= HyperedgeConnectedComponents[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]
Out[1]= {{1, 2, 3}, {4}}

In[2]:= HyperedgeConnectedComponents[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 2]
Out[2]= {{1, 2, 3}, {4}}

In[3]:= HyperedgeConnectedComponents[Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}], 3]
Out[3]= {{1}, {2}, {3}}

In[4]:= HyperedgeConnectedComponents[Hypergraph[{1,2},{{1,2},{},{1}}]]
Out[4]= {{1, 3}}

In[5]:= HyperedgeConnectedComponents[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 0]
Out[5]= HyperedgeConnectedComponents[Hypergraph[<7 vertices, 4 hyperedges>], 0]
```

### Applications (3)

```mathematica
In[6]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
Out[6]= Hypergraph[<7 vertices, 4 hyperedges>]
```

Hyperedge components, as 1-based index lists

```mathematica
In[7]:= HyperedgeConnectedComponents[h]
Out[7]= {{1, 2, 3}, {4}}
```

With s = 2, joined only on a 2-vertex overlap

```mathematica
In[8]:= HyperedgeConnectedComponents[Hypergraph[{{1, 2, 3}, {2, 3, 4}, {5, 6}}], 2]
Out[8]= {{1, 2}, {3}}
```

## Implementation notes

**Algorithm.** `builtin_hyperedge_connected_components` groups the hyperedges
under intersection, returning Lists of 1-based hyperedge indices. The default
`s = 1` runs union–find directly over the vertex→hyperedge incidence lists
(hyperedges sharing a vertex are unioned). For `s >= 2` (the s-connectivity of
Aksoy et al.) each hyperedge is joined to the `s_neighbours` it shares at least
`s` vertices with, then the classes are closed transitively by union–find. A
hyperedge with fewer than `s` distinct vertices lies on no s-walk and is marked
inactive, so it is omitted (empty hyperedges never appear). `groups_to_list`
orders classes by smallest index.

**Data structures.** The distinct-vertex CSR `soff/sv` and incidence CSR
`voff/ve`; a union–find parent array; an `active` mask gating hyperedges by arity;
for `s >= 2` the `cnt`/`touched`/`nb` overlap scratch. The result is a `List` of
index Lists.

**Complexity / limits.** `s = 1` is near-linear in the total incidence; `s >= 2`
is `O(Σ_v deg(v)²)` (overlap counting). `s` must be a positive integer, else the
call is left unevaluated.

- s-connectivity follows Aksoy et al. (2020). The default is `s = 1`.
- Hyperedges with fewer than `s` vertices lie on no s-walk and are omitted (so
  empty hyperedges never appear).
- `s = 1` is union–find over the incidence lists; `s ≥ 2` counts overlaps per
  hyperedge, `O(Σ_v deg(v)²)`.
- `s` must be a positive integer; otherwise unevaluated.
- Benchmark (experiment 96): `s = 2` on 10⁵ hyperedges 18.4 ms warm / 18.6 ms
  cold (Mathematica 1306 ms, xgi 2466 ms).

**Attributes:** `Protected`.

## References

- S. G. Aksoy, C. Joslyn, C. Ortiz Marrero, B. Praggastis and E. Purvine, *Hypernetwork science via high-order hypergraph walks*, EPJ Data Science **9**:16 (2020).
- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

Where `HypergraphConnectedComponents` groups vertices, `HyperedgeConnectedComponents`
groups the **hyperedges**: two hyperedges are in the same component when a chain of
pairwise-intersecting hyperedges links them. Components come back as Lists of
1-based hyperedge indices (`EdgeList` positions), ordered by smallest index.

`HyperedgeConnectedComponents[h, s]` gives the s-connected components of Aksoy et
al.: consecutive hyperedges must share at least `s` vertices. Hyperedges with fewer
than `s` distinct vertices lie on no s-walk and are omitted, so empty hyperedges
never appear. `s = 1` is union–find; `s >= 2` counts overlaps. `s` must be a
positive integer.

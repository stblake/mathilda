# HyperedgeDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HyperedgeDistance[h, i, j] gives the length of a shortest walk of intersecting hyperedges from hyperedge i to hyperedge j, or Infinity. HyperedgeDistance[h, i, j, s] uses s-walks (consecutive hyperedges share at least s vertices). HyperedgeDistance[h, i] or [h, i, All, s] gives the distances from i to every hyperedge.`**

## Examples (6)

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

## Implementation notes

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

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

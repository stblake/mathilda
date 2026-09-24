# HyperedgeConnectedComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HyperedgeConnectedComponents[h] gives the components of the hyperedges of h under intersection, as Lists of 1-based hyperedge indices. HyperedgeConnectedComponents[h, s] gives the s-connected components: hyperedges with at least s vertices, joined when they share at least s.`**

## Examples (5)

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

## Implementation notes

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

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

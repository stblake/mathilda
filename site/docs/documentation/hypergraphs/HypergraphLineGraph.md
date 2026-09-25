# HypergraphLineGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphLineGraph[h] gives the line graph of h: vertices 1..m, i<->j when hyperedges i and j intersect. HypergraphLineGraph[h, s] gives the s-line graph, joining hyperedges that share at least s vertices.`**

## Examples (4)

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

## Implementation notes

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

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

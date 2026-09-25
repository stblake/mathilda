# HypergraphStarExpansion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphStarExpansion[h] gives the incidence (star) expansion of h: the bipartite Graph on VertexList[h] and nodes Hyperedge[1], ..., Hyperedge[m], with v<->Hyperedge[j] whenever v lies in hyperedge j.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= InputForm[HypergraphStarExpansion[Hypergraph[{{1,2},{2,3}}]]]
Out[1]= Graph[{1, 2, 3, Hyperedge[1], Hyperedge[2]}, {1 <-> Hyperedge[1], 2 <-> Hyperedge[1], 2 <-> Hyperedge[2], 3 <-> Hyperedge[2]}]

In[2]:= EdgeCount[HypergraphStarExpansion[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[2]= 9

In[3]:= InputForm[HypergraphStarExpansion[Hypergraph[{Hyperedge[1], 2},{{Hyperedge[1],2}}]]]
Out[3]= HypergraphStarExpansion[Hypergraph[{Hyperedge[1], 2}, {{Hyperedge[1], 2}}]]
```

## Implementation notes

- Vertices are `VertexList[h]` followed by nodes `Hyperedge[1], ...,
  Hyperedge[m]`, with an edge `v <-> Hyperedge[j]` for each `v ∈ e_j`.
- Returns a simple `Graph`, validated and memoized as usual.
- Unevaluated if some vertex of `h` is itself such a `Hyperedge[j]` node.
- Benchmark (experiment 96): 10⁵ hyperedges 38.5 ms warm / 38.3 ms cold
  (Mathematica 313 ms, xgi 466 ms).

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/)

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

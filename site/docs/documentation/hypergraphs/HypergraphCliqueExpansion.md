# HypergraphCliqueExpansion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphCliqueExpansion[h] gives the 2-section of h: the undirected Graph on VertexList[h] with u<->v whenever u and v share a hyperedge.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= InputForm[HypergraphCliqueExpansion[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[1]= Graph[{1, 2, 3, 4, 5, 6, 7}, {1 <-> 2, 1 <-> 3, 2 <-> 3, 3 <-> 4, 4 <-> 5, 4 <-> 6, 5 <-> 6}]

In[2]:= InputForm[HypergraphCliqueExpansion[Hypergraph[{{1,1,2},{2,3}}]]]
Out[2]= Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}]

In[3]:= HypergraphCliqueExpansion[{{1, 2}}]
Out[3]= Graph[<2 vertices, 1 edge>]
```

## Implementation notes

- Returns a simple `Graph` (validated and memoized as usual) on `VertexList[h]`;
  edges in order of first co-occurrence. A vertex repeated inside a hyperedge
  gives no self-loop.
- Also equals the graph whose distance `HypergraphDistance` measures.
- Benchmark (experiment 96): 10⁵ hyperedges 32.3 ms warm / 33.8 ms cold
  (Mathematica 278 ms, xgi 729 ms).

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/), [HypergraphDistance](../../hypergraphs/HypergraphDistance/)

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

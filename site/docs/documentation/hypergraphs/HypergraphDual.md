# HypergraphDual

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphDual[h] gives the dual hypergraph: vertices 1..m (one per hyperedge of h) and, for each vertex v of h in VertexList order, the hyperedge of indices of the hyperedges containing v.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= InputForm[HypergraphDual[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[1]= Hypergraph[{1, 2, 3, 4}, {{1}, {1}, {1, 2}, {2, 3}, {3}, {3}, {4}}]

In[2]:= InputForm[HypergraphDual[HypergraphDual[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]]
Out[2]= Hypergraph[{1, 2, 3, 4, 5, 6, 7}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]

In[3]:= InputForm[HypergraphDual[Hypergraph[{a,b,c},{{a,b},{a,b}}]]]
Out[3]= Hypergraph[{1, 2}, {{1, 2}, {1, 2}, {}}]

In[4]:= HypergraphDual[{{1, 2}}]
Out[4]= HypergraphDual[{{1, 2}}]
```

## Implementation notes

- Vertices are `1..m`; for each vertex of `h`, in VertexList order, the dual has
  the hyperedge of the (ascending) indices of the hyperedges containing it. An
  isolated vertex gives an empty hyperedge.
- Hyperedges are read as sets. The dual of the dual recovers the incidence
  structure on `1..n`.
- Linear in the total incidence. Benchmark (experiment 96): 10⁵ hyperedges
  8.1 ms warm / 10.1 ms cold (Mathematica 236 ms, xgi 361 ms).
- A bare List of hyperedges is not accepted (unevaluated).

**Attributes:** `Protected`.

## References

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

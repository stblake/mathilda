# ConnectedHypergraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ConnectedHypergraphQ[h] gives True if h has at least one vertex and is connected. h may be a plain List of hyperedges.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= HypergraphConnectedComponents[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]
Out[1]= {{1, 2, 3, 4, 5, 6}, {7}}

In[2]:= HypergraphConnectedComponents[Hypergraph[{c,b,a},{{a,c}}]]
Out[2]= {{c, a}, {b}}

In[3]:= ConnectedHypergraphQ[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]
Out[3]= False

In[4]:= ConnectedHypergraphQ[{{1,2},{2,3}}]
Out[4]= True

In[5]:= ConnectedHypergraphQ[{}]
Out[5]= False
```

## Implementation notes

- Components are ordered by their first vertex (Mathilda's
  `ConnectedComponents` convention), vertices within a component in VertexList
  order. Isolated vertices are singleton components.
- Union–find, near-linear in the total incidence.
- `ConnectedHypergraphQ` is the Wolfram Function Repository name. It accepts a
  bare List of hyperedges; the FR function leaves `{}` unevaluated, Mathilda
  gives `False`. It gives `False` for any non-hypergraph.
- Benchmark (experiment 96): `HypergraphConnectedComponents` on 10⁵ hyperedges
  2.3 ms warm and cold (Mathematica 160 ms, xgi 282 ms);
  `ConnectedHypergraphQ` 20.1 ms (FR function 261 ms, xgi 546 ms).

**Attributes:** `Protected`.

## References

**See also:** [HypergraphConnectedComponents](../../hypergraphs/HypergraphConnectedComponents/), [ConnectedComponents](../../graphs/ConnectedComponents/)

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

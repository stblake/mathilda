# HypergraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphQ[h] gives True if h is a valid Hypergraph, and False otherwise.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= HypergraphQ[Hypergraph[{{1,2,3},{3,4}}]]
Out[1]= True

In[2]:= HypergraphQ[Graph[{1<->2}]]
Out[2]= False

In[3]:= GraphQ[Hypergraph[{{1,2,3},{3,4}}]]
Out[3]= False

In[4]:= HypergraphQ[Hypergraph[{1,2},{{1,5}}]]
Out[4]= False
```

## Implementation notes

- `False` for a `Graph`, for a bare List of hyperedges, and for a malformed
  (unevaluated) `Hypergraph[...]`. Conversely `GraphQ` of a hypergraph is `False`.
- `O(1)` on a memoized hypergraph object.

**Attributes:** `Protected`.

## References

**See also:** [Hypergraph](../../hypergraphs/Hypergraph/), [Graph](../../graphs/Graph/), [GraphQ](../../graphs/GraphQ/)

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

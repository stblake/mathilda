# HypergraphEdgeDelete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphEdgeDelete[h, e] or [h, {e1, ...}] deletes every hyperedge identical (SameQ) to one given.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= InputForm[HypergraphEdgeAdd[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {7, 8}]]
Out[1]= Hypergraph[{1, 2, 3, 4, 5, 6, 7, 8}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}, {7, 8}}]

In[2]:= InputForm[HypergraphEdgeAdd[Hypergraph[{{1,2}}], {1,2}]]
Out[2]= Hypergraph[{1, 2}, {{1, 2}, {1, 2}}]

In[3]:= InputForm[HypergraphEdgeDelete[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {{3,4},{7}}]]
Out[3]= Hypergraph[{1, 2, 3, 4, 5, 6, 7}, {{1, 2, 3}, {4, 5, 6}}]

In[4]:= InputForm[HypergraphEdgeDelete[Hypergraph[{{1,2},{2,3},{1,2}}], {1,2}]]
Out[4]= Hypergraph[{1, 2, 3}, {{2, 3}}]

In[5]:= HypergraphEdgeDelete[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {4, 3}]
Out[5]= HypergraphEdgeDelete[Hypergraph[<7 vertices, 4 hyperedges>], {4, 3}]
```

## Implementation notes

- In the Edge heads, a List whose every element is a List is a list of
  hyperedges; otherwise it is one hyperedge.
- `HypergraphEdgeAdd` allows repeats (a multi-hypergraph).
- `HypergraphEdgeDelete` compares hyperedges as written (`{2, 1}` does not
  delete `{1, 2}`) and removes all copies of a repeated hyperedge. It is
  unevaluated if a named hyperedge does not occur.

**Attributes:** `Protected`.

## References

**See also:** [HypergraphEdgeAdd](../../hypergraphs/HypergraphEdgeAdd/), [SameQ](../../comparisons/SameQ/)

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

# HypergraphToGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphToGraph[h] converts h, read as an ordered hypergraph, to the directed Graph with v_a->v_b for every a<b in each hyperedge {v_1, ..., v_k} (the Function Repository's HypergraphToGraph). Self-loops are dropped and parallel edges merged. h may be a plain List of hyperedges.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= InputForm[HypergraphToGraph[{{1,2,3},{3,4}}]]
Out[1]= Graph[{1, 2, 3, 4}, {1 -> 2, 1 -> 3, 2 -> 3, 3 -> 4}]

In[2]:= InputForm[HypergraphToGraph[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[2]= Graph[{1, 2, 3, 4, 5, 6, 7}, {1 -> 2, 1 -> 3, 2 -> 3, 3 -> 4, 4 -> 5, 4 -> 6, 5 -> 6}]

In[3]:= InputForm[HypergraphToGraph[{{1,1,2},{1,2},{5}}]]
Out[3]= Graph[{1, 2, 5}, {1 -> 2}]

In[4]:= HypergraphToGraph[5]
Out[4]= HypergraphToGraph[5]
```

## Implementation notes

- Wolfram Function Repository name and semantics: hyperedge `{v1, ..., vk}`
  contributes `v_a -> v_b` for every `a < b`.
- Accepts a bare List of hyperedges, as the FR function does.
- *Deviation:* Mathilda graphs are simple, so self-loops (from a repeated
  vertex) are dropped and parallel copies merged, and every vertex is kept (the
  FR function returns a multigraph and drops vertices that occur only in unary
  hyperedges).
- Benchmark (experiment 96): 10⁵ hyperedges 59.3 ms (FR function 233 ms, xgi
  376 ms).

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/)

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

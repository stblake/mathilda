# Subhypergraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Subhypergraph[h, {v1, ...}] gives the sub-hypergraph induced by the given vertices: those of them in h, and the hyperedges lying entirely among them.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= InputForm[Subhypergraph[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {1, 2, 3, 4}]]
Out[1]= Hypergraph[{1, 2, 3, 4}, {{1, 2, 3}, {3, 4}}]

In[2]:= InputForm[Subhypergraph[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {4, 3, 99}]]
Out[2]= Hypergraph[{3, 4}, {{3, 4}}]

In[3]:= InputForm[HypergraphRestriction[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {1, 2, 3, 4}]]
Out[3]= Hypergraph[{1, 2, 3, 4}, {{1, 2, 3}, {3, 4}, {4}}]

In[4]:= InputForm[HypergraphRestriction[Hypergraph[{{3,1,3,2},{5}}], {3, 2}]]
Out[4]= Hypergraph[{3, 2}, {{3, 3, 2}}]

In[5]:= Subhypergraph[{{1,2}}, {1}]
Out[5]= Subhypergraph[{{1, 2}}, {1}]
```

## Implementation notes

- `Subhypergraph` follows `Subgraph`'s rule: the named vertices that are in `h`
  (in VertexList order; absent names are ignored) and the hyperedges lying
  entirely among them.
- `HypergraphRestriction` preserves the order and repeats of the kept vertices
  within each hyperedge, and drops hyperedges that miss the vertex set.
- A bare List of hyperedges is not accepted (unevaluated).

**Attributes:** `Protected`.

## References

**See also:** [HypergraphRestriction](../../hypergraphs/HypergraphRestriction/), [Subgraph](../../graphs/Subgraph/)

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

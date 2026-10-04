# Subhypergraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Subhypergraph[h, {v1, ...}] gives the sub-hypergraph induced by the given vertices: those of them in h, and the hyperedges lying entirely among them.`**

## Examples (8)

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

### Applications (3)

```mathematica
In[6]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
Out[6]= Hypergraph[<7 vertices, 4 hyperedges>]
```

Keep only hyperedges lying entirely inside

```mathematica
In[7]:= InputForm[Subhypergraph[h, {1, 2, 3, 4}]]
Out[7]= Hypergraph[{1, 2, 3, 4}, {{1, 2, 3}, {3, 4}}]
```

{3, 4} straddles the boundary and is dropped

```mathematica
In[8]:= InputForm[Subhypergraph[h, {4, 5, 6}]]
Out[8]= Hypergraph[{4, 5, 6}, {{4, 5, 6}}]
```

## Implementation notes

**Algorithm.** `builtin_subhypergraph` follows `Subgraph`'s rule: keep the named
vertices that are in `h` (in `VertexList` order; absent names ignored) and the
hyperedges lying **entirely** among them. `vertex_mask` turns the vertex-List
argument into a membership bitmap `in` over `h`'s vertices (non-vertices ignored).
A hyperedge is kept only if every vertex of its distinct set is in the mask; then
`rebuild` emits the kept vertices and hyperedges in original order.

**Data structures.** The distinct-vertex CSR `soff/sv`, the memo's vertex
`GraphVIdx` (for `vertex_mask`), an `in` bitmap and an `ke` hyperedge keep-mask;
`rebuild` constructs the fresh `Hypergraph`.

**Complexity / limits.** `O(n + Σ|e|)`. Contrast `HypergraphRestriction`, which
keeps *partial* hyperedges by intersecting rather than requiring full containment.
A bare List of hyperedges is not accepted (left unevaluated).

- `Subhypergraph` follows `Subgraph`'s rule: the named vertices that are in `h`
  (in VertexList order; absent names are ignored) and the hyperedges lying
  entirely among them.
- `HypergraphRestriction` preserves the order and repeats of the kept vertices
  within each hyperedge, and drops hyperedges that miss the vertex set.
- A bare List of hyperedges is not accepted (unevaluated).

**Attributes:** `Protected`.

## References

**See also:** [HypergraphRestriction](../../hypergraphs/HypergraphRestriction/), [Subgraph](../../graphs/Subgraph/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`Subhypergraph[h, {v1, ...}]` follows `Subgraph`'s rule: it keeps the named
vertices that occur in `h` (in `VertexList` order; absent names are ignored) and
only the hyperedges lying **entirely** among them. A hyperedge with even one
vertex outside the set is dropped whole.

Contrast `HypergraphRestriction`, which keeps such a hyperedge in truncated form by
intersecting it with the vertex set. A bare List of hyperedges is not accepted
here. The operation is linear in the total incidence.

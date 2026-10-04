# HypergraphRestriction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphRestriction[h, {v1, ...}] intersects every hyperedge of h with the given vertices, dropping hyperedges that miss them (Berge's induced sub-hypergraph).`**

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

Each hyperedge intersected with {1,2,3,4}

```mathematica
In[7]:= InputForm[HypergraphRestriction[h, {1, 2, 3, 4}]]
Out[7]= Hypergraph[{1, 2, 3, 4}, {{1, 2, 3}, {3, 4}, {4}}]
```

Hyperedges missing the set entirely are dropped

```mathematica
In[8]:= InputForm[HypergraphRestriction[h, {1, 2}]]
Out[8]= Hypergraph[{1, 2}, {{1, 2}}]
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_restriction` gives Berge's induced
sub-hypergraph: every hyperedge of `h` intersected with the given vertex set, with
hyperedges that miss it entirely dropped. `vertex_mask` builds the membership
bitmap `in`. Each hyperedge is walked over its **raw** elements `ev` (so order and
repeats of the kept vertices survive); if all its elements are kept the original
`List` is shared unchanged, otherwise a new `List` of the kept elements is built.
An intersection that is empty is dropped. The kept vertices (in `VertexList` order)
form the new vertex set.

**Data structures.** The raw hyperedge CSR `eoff/ev` and distinct set are read
through the `HypView`; an `in` bitmap, and `vs`/`es`/`tmp` buffers for the rebuilt
vertices, hyperedges, and per-hyperedge kept elements. The result is a fresh
`Hypergraph`.

**Complexity / limits.** `O(n + Σ|e|)`. Unlike `Subhypergraph`, a hyperedge
straddling the boundary is kept in truncated form, so the result may have more
hyperedges than `Subhypergraph` on the same vertex set. A bare List is not
accepted.

- `Subhypergraph` follows `Subgraph`'s rule: the named vertices that are in `h`
  (in VertexList order; absent names are ignored) and the hyperedges lying
  entirely among them.
- `HypergraphRestriction` preserves the order and repeats of the kept vertices
  within each hyperedge, and drops hyperedges that miss the vertex set.
- A bare List of hyperedges is not accepted (unevaluated).

**Attributes:** `Protected`.

## References

**See also:** [Subhypergraph](../../hypergraphs/Subhypergraph/), [Subgraph](../../graphs/Subgraph/)

- C. Berge, *Hypergraphs: Combinatorics of Finite Sets*, North-Holland Mathematical Library 45 (Elsevier, 1989), ch. 1 (induced sub-hypergraph).
- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`HypergraphRestriction[h, {v1, ...}]` is Berge's induced sub-hypergraph: every
hyperedge of `h` is intersected with the given vertex set, and hyperedges that miss
it entirely are dropped. Restricting the running example to `{1, 2, 3, 4}` keeps
`{1, 2, 3}` and `{3, 4}` whole but shrinks `{4, 5, 6}` to `{4}`.

The order and repeats of the kept vertices are preserved within each hyperedge.
This is the difference from `Subhypergraph`, which instead keeps only the
hyperedges that fit entirely inside the set, discarding the rest; restriction can
therefore return more hyperedges than `Subhypergraph` on the same vertices. A bare
List of hyperedges is not accepted.

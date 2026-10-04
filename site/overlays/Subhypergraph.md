### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= InputForm[Subhypergraph[h, {1, 2, 3, 4}]]  (* keep only hyperedges lying entirely inside *)
```

```mathematica
In[1]:= InputForm[Subhypergraph[h, {4, 5, 6}]]  (* {3, 4} straddles the boundary and is dropped *)
```

### Notes

`Subhypergraph[h, {v1, ...}]` follows `Subgraph`'s rule: it keeps the named
vertices that occur in `h` (in `VertexList` order; absent names are ignored) and
only the hyperedges lying **entirely** among them. A hyperedge with even one
vertex outside the set is dropped whole.

Contrast `HypergraphRestriction`, which keeps such a hyperedge in truncated form by
intersecting it with the vertex set. A bare List of hyperedges is not accepted
here. The operation is linear in the total incidence.

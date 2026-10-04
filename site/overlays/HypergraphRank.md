### Worked examples

```mathematica
In[1]:= HypergraphRank[Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]]  (* the largest arity *)
```

```mathematica
In[1]:= HypergraphRank[Hypergraph[{1, 2}, {}]]  (* no hyperedges -> 0 *)
```

### Notes

`HypergraphRank[h]` is the largest hyperedge arity — `Max[HyperedgeSizes[h]]` — and
`HypergraphCorank[h]` the smallest. Arity is the `Length` as written, so a repeated
vertex contributes to the rank. A hypergraph with no hyperedges has rank `0`.

A hypergraph is `k`-uniform exactly when its rank equals its corank equals `k`;
`UniformHypergraphQ` tests that directly. Like the other arity heads this reads
only the hyperedge offsets, so it is `O(m)` and independent of how large the
hyperedges are.

### Worked examples

```mathematica
In[1]:= HypergraphCorank[Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]]  (* the smallest arity *)
```

```mathematica
In[1]:= HypergraphCorank[Hypergraph[{{1, 2, 3}, {2, 3, 4}}]]  (* a 3-uniform hypergraph: corank = rank *)
```

### Notes

`HypergraphCorank[h]` is the smallest hyperedge arity — `Min[HyperedgeSizes[h]]` —
the companion of `HypergraphRank` (the largest). Arity is the `Length` as written,
counting a repeated vertex; with no hyperedges the corank is `0`.

When rank and corank coincide the hypergraph is uniform, and their common value is
the `k` of `k`-uniformity. The computation reads only the hyperedge offsets, so it
is `O(m)`.

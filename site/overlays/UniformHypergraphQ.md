### Worked examples

```mathematica
In[1]:= UniformHypergraphQ[Hypergraph[{{1, 2, 3}, {2, 3, 4}}]]  (* all hyperedges have arity 3 *)
```

```mathematica
In[1]:= UniformHypergraphQ[Hypergraph[{{1, 2, 3}, {3, 4}}]]  (* mixed arities *)
```

```mathematica
In[1]:= UniformHypergraphQ[Hypergraph[{{1, 2, 3}, {2, 3, 4}}], 2]  (* 3-uniform, so not 2-uniform *)
```

### Notes

`UniformHypergraphQ[h]` is `True` when every hyperedge has the same arity;
`UniformHypergraphQ[h, k]` additionally pins that common arity to `k`. Arity is the
`Length` as written, so repeats count toward uniformity, and a hypergraph with no
hyperedges is vacuously uniform.

Equivalently, `h` is uniform when `HypergraphRank[h] == HypergraphCorank[h]`. As a
`*Q` predicate the one-argument form returns `False` on a non-hypergraph rather
than staying unevaluated; the `k` form with a malformed `k` is left unevaluated.

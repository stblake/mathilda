### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= InputForm[HypergraphRestriction[h, {1, 2, 3, 4}]]  (* each hyperedge intersected with {1,2,3,4} *)
```

```mathematica
In[1]:= InputForm[HypergraphRestriction[h, {1, 2}]]  (* hyperedges missing the set entirely are dropped *)
```

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

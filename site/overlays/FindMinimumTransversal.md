### Worked examples

```mathematica
In[1]:= FindMinimumTransversal[{{1, 2}, {2, 3}}]  (* one vertex, 2, hits both hyperedges *)
```

```mathematica
In[1]:= FindMinimumTransversal[{{1, 2}, {2, 3}, {3, 4}}]  (* a smallest hitting set has size 2 *)
```

```mathematica
In[1]:= Length[FindMinimumTransversal[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]]  (* the minimum transversal number *)
```

### Notes

`FindMinimumTransversal[h]` returns a single smallest vertex set meeting every
hyperedge — a minimum hitting set, so its `Length` is the transversal number of
`h`. Where `TransversalHypergraph` enumerates all the *minimal* transversals,
this returns just one of *minimum* cardinality.

The minimum hitting set is NP-hard, so the engine is exact branch and bound: a
greedy max-coverage upper bound seeds the search, which is then pruned by a degree
bound, a disjoint-hyperedge packing bound, and a fractional LP-dual bound — and
when a packing already matches the greedy bound the answer is proved optimal with
no search. Like `TransversalHypergraph` it is correct-or-unevaluated: it is left
unevaluated if some hyperedge is empty (nothing can hit it) or if the exact search
exceeds its node budget, and `TimeConstrained` can interrupt it. Both a
`Hypergraph` and a bare List of hyperedges are accepted.

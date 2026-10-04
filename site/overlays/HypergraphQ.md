### Worked examples

```mathematica
In[1]:= HypergraphQ[Hypergraph[{{1, 1, 2}, {2, 3}}]]  (* a repeated vertex inside a hyperedge is allowed *)
```

```mathematica
In[1]:= GraphQ[Hypergraph[{{1, 2}}]]  (* a hypergraph is not a Graph *)
```

```mathematica
In[1]:= HypergraphQ["not a hypergraph"]
```

### Notes

`HypergraphQ` is `True` only for a canonical, valid `Hypergraph` object: a
`Hypergraph[List, List]` with pairwise-distinct vertices and every hyperedge a
`List` of declared vertices. It is `False` for a `Graph`, for a bare List of
hyperedges, and for a malformed `Hypergraph[...]` that the constructor left
unevaluated (for instance one naming a vertex that is not in its vertex List).

The check is `O(1)` once the object has been seen, because validation is the same
pointer-keyed memo lookup every hypergraph head shares. Like the other `*Q`
predicates it always returns a Boolean and never leaves itself unevaluated.

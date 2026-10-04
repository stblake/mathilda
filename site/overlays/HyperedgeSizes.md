### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= HyperedgeSizes[h]  (* the arity of each hyperedge, in EdgeList order *)
```

```mathematica
In[1]:= Total[HyperedgeSizes[h]]  (* their total is the total incidence of h *)
```

```mathematica
In[1]:= HyperedgeSizes[Hypergraph[{{1, 1, 2}, {}}]]  (* arity counts repeats; an empty hyperedge is 0 *)
```

### Notes

`HyperedgeSizes` reports each hyperedge's `Length` **as written** — a repeated
vertex counts, and an empty hyperedge is `0` — because arity is an ordered,
multiset notion, not the distinct-vertex count. The sum of the sizes is the total
incidence `Σ|e|`, which is also the edge count of the star expansion.

The result is a packed `int64` list above the packing threshold, so it is cheap to
feed into `Max`, `Min`, `Total`, `Tally` and the rest of the numeric pipeline —
indeed `HypergraphRank` and `HypergraphCorank` are just its maximum and minimum.

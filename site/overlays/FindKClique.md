### Worked examples

```mathematica
In[1]:= FindKClique[Graph[{1 <-> 2, 2 <-> 3, 1 <-> 3, 3 <-> 4}], 1]  (* k = 1 is an ordinary maximum clique *)
```

```mathematica
In[1]:= FindKClique[PathGraph[5], 2]  (* k = 2: vertices pairwise within distance 2 *)
```

### Notes

`FindKClique[g, k]` returns `{c}`, where `c` is a largest set of vertices that are pairwise
within graph distance `k`. With `k = 1` this is a maximum clique of `g`; larger `k` relaxes
adjacency to "close enough", so on a path the three consecutive vertices within distance 2 of
each other form the answer.

The empty graph gives `{}`. The search is exact: if the node budget is exhausted the call stays
unevaluated rather than return a set that might not be largest.

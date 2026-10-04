### Worked examples

```mathematica
In[1]:= BetweennessCentrality[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5}]]  (* the interior of a path carries every cross-path *)
```

```mathematica
In[1]:= BetweennessCentrality[Graph[{1 -> 2, 2 -> 3, 1 -> 3, 3 -> 4}]]  (* directed: ordered pairs, no 1/2 factor *)
```

### Notes

The value is unnormalised Freeman betweenness: for each vertex `v`, the sum over pairs of
other vertices `s`, `t` of the fraction of shortest `s`-`t` paths that pass through `v`.
Undirected graphs count each pair once; directed graphs count ordered pairs.

The vertex form ignores `EdgeWeight` and always measures shortest paths by edge count. If you
need weighted shortest paths for the pair fractions, that distinction lives in
`EdgeBetweennessCentrality`, which does honour `EdgeWeight`. A mixed graph (both directed and
undirected edges) is left unevaluated.

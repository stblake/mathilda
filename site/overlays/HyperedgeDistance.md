### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= HyperedgeDistance[h, 1, 3]  (* shortest walk of intersecting hyperedges, 1 -> 2 -> 3 *)
```

```mathematica
In[1]:= HyperedgeDistance[h, 1]  (* distances from hyperedge 1 to every hyperedge *)
```

```mathematica
In[1]:= HyperedgeDistance[h, 1, 4]  (* hyperedge 4 = {7} is unreachable *)
```

### Notes

`HyperedgeDistance[h, i, j]` is the number of steps in a shortest walk of
intersecting hyperedges from `i` to `j`, `0` when `i == j` and `Infinity` when no
walk exists. `HyperedgeDistance[h, i, j, s]` uses s-walks (consecutive hyperedges
sharing at least `s` vertices); `HyperedgeDistance[h, i]` and
`HyperedgeDistance[h, i, All, s]` give the distances from `i` to every hyperedge.

The search is a BFS run directly on the incidence structure — the line graph is
never built — so the single-source form is linear in the total incidence for
`s = 1`. Hyperedge indices are 1-based `EdgeList` positions; an out-of-range index
leaves the call unevaluated.

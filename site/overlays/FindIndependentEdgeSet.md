### Worked examples

```mathematica
In[1]:= FindIndependentEdgeSet[CompleteGraph[4]]  (* the complete graph on 4 vertices has a perfect matching of two edges *)
```

```mathematica
In[1]:= FindIndependentEdgeSet[PathGraph[{1, 2, 3, 4}]]  (* alternate edges of the path *)
```

```mathematica
In[1]:= FindIndependentEdgeSet[CycleGraph[5]]  (* an odd cycle matches all but one vertex *)
```

### Notes

`FindIndependentEdgeSet[g]` returns a maximum **matching** — a largest set of
pairwise non-adjacent edges — as `g`'s own edge expressions in `EdgeList` order.
Edge direction is ignored. The size of the result is the matching number
`ν(g)`; the particular maximum matching returned is not specified.

The solver is exact: a Karp-Sipser greedy start, then Hopcroft-Karp for
bipartite graphs and Edmonds' blossom algorithm for general ones. Because K4 and
the 4-path both admit a perfect (resp. near-perfect) matching, the odd 5-cycle
can cover only four of its five vertices, leaving one exposed.

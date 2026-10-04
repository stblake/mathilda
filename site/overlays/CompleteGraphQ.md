### Worked examples

```mathematica
In[1]:= CompleteGraphQ[CompleteGraph[4]]  (* the complete graph on 4 vertices *)
```

```mathematica
In[1]:= CompleteGraphQ[CycleGraph[5]]  (* a 5-cycle is far from complete *)
```

```mathematica
In[1]:= CompleteGraphQ[CompleteGraph[4], {1, 2, 3}]  (* the subgraph induced by three vertices is still complete *)
```

```mathematica
In[1]:= CompleteGraphQ[CycleGraph[5], {1, 2, 3}]  (* those three vertices are not pairwise adjacent *)
```

```mathematica
In[1]:= CompleteGraphQ[5]  (* a non-graph argument is simply False *)
```

### Notes

Completeness means every pair of distinct vertices is joined; for a directed graph that
requires edges in *both* directions (equivalently, an undirected edge). Graphs on zero or
one vertices are complete by default, since there is no pair that could fail.

The two-argument form `CompleteGraphQ[g, vlist]` tests the subgraph of `g` induced by
`vlist`, and gives `False` if any element of `vlist` is not actually a vertex of `g`.
Like every `*Q` predicate, a non-graph argument yields `False` rather than staying
unevaluated.

### Worked examples

```mathematica
In[1]:= FindVertexColoring[CompleteGraph[4]]  (* a clique needs a distinct colour per vertex *)
```

```mathematica
In[1]:= FindVertexColoring[CycleGraph[4]]  (* an even cycle is 2-chromatic *)
```

```mathematica
In[1]:= FindVertexColoring[CycleGraph[5]]  (* an odd cycle needs three colours *)
```

```mathematica
In[1]:= FindVertexColoring[PathGraph[{1, 2, 3, 4}]]  (* a path is bipartite *)
```

```mathematica
In[1]:= Max[FindVertexColoring[PetersenGraph[]]]  (* the chromatic number of the Petersen graph is 3 *)
```

### Notes

`FindVertexColoring[g]` returns a list of positive-integer colour labels, one per
vertex in `VertexList` order, such that adjacent vertices differ — and, crucially,
using as few distinct colours as possible. The number of distinct colours is the
chromatic number of `g`, so `Max` of the result reads off that number. An edge
constrains its endpoints in either direction, so direction is ignored.

Because minimal colouring is NP-hard, the search is exact and may **refuse** rather
than return a merely-valid colouring: the call is left unevaluated for a graph of
more than 128 vertices, or when an internal node budget is exhausted before
minimality can be proven. Bound an interactive call with `TimeConstrained` if
needed. Only the one-argument form is supported.

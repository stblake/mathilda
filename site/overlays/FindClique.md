### Worked examples

```mathematica
In[1]:= FindClique[CompleteGraph[4]]  (* the whole graph is one clique *)
```

```mathematica
In[1]:= FindClique[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 1, 3 <-> 4}]]  (* the triangle, not the pendant edge *)
```

```mathematica
In[1]:= FindClique[CycleGraph[5]]  (* triangle-free: the maximum clique is a single edge *)
```

```mathematica
In[1]:= FindClique[CompleteGraph[4], 2]  (* maximal cliques of <= 2 vertices -- none, every maximal clique is K4 *)
```

```mathematica
In[1]:= FindClique[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 1, 3 <-> 4}], {2}, All]  (* every maximal clique of exactly 2 vertices *)
```

### Notes

`FindClique[g]` returns a one-element list `{c}` holding one maximum clique, each
clique in `VertexList` order. The maximum is proved by an exact bitset
branch-and-bound (Tomita MCS with a bit-parallel colouring bound); which maximum
clique is returned is not specified, as in Mathematica.

The size-spec forms report **maximal** cliques — not extendable — so
`FindClique[g, 2]` is empty whenever every maximal clique is larger than two
vertices. A spec may be `k` (at most `k`), `{k}` (exactly `k`) or `{kmin, kmax}`,
and a trailing count or `All` asks for several, largest first. In a directed
graph a clique needs the arcs in both directions (an `UndirectedEdge` counts as
both). The search is budgeted and leaves the call unevaluated on exhaustion
rather than returning a non-maximum clique.

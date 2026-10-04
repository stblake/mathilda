### Worked examples

```mathematica
In[1]:= g1 = Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4}];
```

```mathematica
In[1]:= g2 = Graph[{2 <-> 3, 3 <-> 4, 4 <-> 1}];
```

```mathematica
In[1]:= EdgeList[GraphIntersection[g1, g2]]  (* the edges shared by both *)
```

### Notes

`GraphIntersection[g1, g2, ...]` gives the graph on the union of the vertex sets whose edges are
those common to all of the `gi`, in canonical order. Weights are dropped.

An undirected edge matches its reversal, so `1 <-> 2` and `2 <-> 1` are the same edge for the
set operation. The result is a canonical `Graph`; query it with `EdgeList` / `EdgeCount`.

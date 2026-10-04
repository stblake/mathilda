### Worked examples

```mathematica
In[1]:= g1 = Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4}];
```

```mathematica
In[1]:= g2 = Graph[{2 <-> 3, 3 <-> 4, 4 <-> 1}];
```

```mathematica
In[1]:= EdgeList[GraphDifference[g1, g2]]  (* edges of g1 that g2 does not have *)
```

### Notes

`GraphDifference[g1, g2]` gives the graph on the union of the vertex sets whose edges are those
of `g1` not present in `g2`, in canonical order. Weights are dropped, and the vertex set is the
union of both graphs' vertices — not just `g1`'s.

An undirected edge matches its reversal for the comparison. The result is a canonical `Graph`;
query it with `EdgeList` / `EdgeCount`.

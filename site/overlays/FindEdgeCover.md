### Worked examples

```mathematica
In[1]:= FindEdgeCover[PathGraph[4]]  (* two edges cover all four vertices *)
```

```mathematica
In[1]:= FindEdgeCover[Graph[{1, 2, 3}, {2 <-> 3}]]  (* an isolated vertex has no edge cover *)
```

### Notes

A minimum edge cover is a smallest set of edges that touches every vertex; its size is `n` minus
the size of a maximum matching. The edges are returned in `EdgeList` order and the set is
guaranteed minimum.

A graph with any isolated vertex has no edge cover at all (nothing can touch that vertex), so the
result is `{}`.

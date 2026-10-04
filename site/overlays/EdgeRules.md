### Worked examples

```mathematica
In[1]:= EdgeRules[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]  (* edges as u -> v rules *)
```

### Notes

`EdgeRules[g]` gives the edges of `g` as a list of rules `u -> v`, in `EdgeList` order. Both
directed and undirected edges become plain `->` rules — the direction is not preserved in the
output, since the point of this form is to feed the edges back into `Graph` or `ReplaceAll`.

The result is an ordinary list of `Rule`s, not a graph; use `EdgeList` if you need the edges with
their `DirectedEdge` / `UndirectedEdge` heads intact.

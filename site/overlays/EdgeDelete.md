### Worked examples

```mathematica
In[1]:= EdgeList[EdgeDelete[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4}], 2 <-> 3]]  (* drop the middle edge *)
```

### Notes

`EdgeDelete[g, e]` removes the edge `e`, keeping every vertex (the vertex set is unchanged).
`EdgeDelete[g, {e1, ...}]` removes several, and `EdgeDelete[g, patt]` removes every edge matching
a pattern. An undirected edge matches either orientation.

A literal edge that is not in `g` leaves the call unevaluated. The result is a canonical `Graph`
object — inspect it with `EdgeList` or `EdgeCount`.

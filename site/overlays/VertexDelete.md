### Worked examples

```mathematica
In[1]:= vd = VertexDelete[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}], 3];
```

```mathematica
In[1]:= {VertexList[vd], EdgeList[vd]}  (* vertex 3 and its three incident edges are gone *)
```

### Notes

`VertexDelete[g, v]` removes vertex `v` and every edge touching it. `VertexDelete[g, {v1, ...}]`
removes several, and `VertexDelete[g, patt]` removes every vertex matching a pattern; the
original vertex and edge order is preserved among the survivors.

Each listed vertex must be a vertex of `g` (or the argument must be a pattern), otherwise the call
stays unevaluated. The result is a canonical `Graph` — inspect it with `VertexList` / `EdgeList`.

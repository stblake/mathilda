### Worked examples

```mathematica
In[1]:= EdgeList[PathGraph[5]]  (* a chain 1 - 2 - 3 - 4 - 5 *)
```

```mathematica
In[1]:= VertexList[PathGraph[{a, b, c}]]  (* a path over the given vertices *)
```

```mathematica
In[1]:= EdgeList[PathGraph[{a, b, c}]]
```

### Notes

`PathGraph[n]` is the undirected path on vertices `1..n`, with `n-1` edges joining consecutive
vertices. `PathGraph[{v1, ..., vk}]` builds the same chain over the given vertex expressions,
which may be any symbols or values.

The result is a canonical `Graph` object; query it with `VertexList`, `EdgeList`, `VertexCount`
or `EdgeCount` rather than relying on its printed form.

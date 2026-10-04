### Worked examples

```mathematica
In[1]:= VertexList[VertexReplace[PathGraph[{1, 2, 3}], 2 -> x]]  (* rename one vertex *)
```

```mathematica
In[1]:= EdgeList[VertexReplace[CycleGraph[4], {1 -> a, 2 -> b}]]  (* edges follow the renamed vertices *)
```

```mathematica
In[1]:= VertexList[VertexReplace[PathGraph[{1, 2, 3}], n_Integer :> n^2]]  (* patterns and delayed rules are allowed *)
```

```mathematica
In[1]:= EdgeList[VertexReplace[Graph[{1 -> 2, 2 -> 3}], {1 -> a, 3 -> c}]]  (* directions are preserved *)
```

```mathematica
In[1]:= EdgeList[VertexReplace[Graph[{1 -> 2, 3 -> 2}], {1 -> 3, 3 -> 1}]]  (* a swap permutes the labels *)
```

### Notes

The rules are applied to the vertex list with the semantics of `Replace` at level 1, so the first matching rule wins and vertices with no match keep their name.

If the renaming sends two vertices to one name and that would create a self-loop or parallel edge, the result would not be a simple graph and the call is left unevaluated.

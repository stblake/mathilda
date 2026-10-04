### Worked examples

```mathematica
In[1]:= VertexQ[CompleteGraph[4], 2]  (* 2 is one of the vertices 1..4 *)
```

```mathematica
In[1]:= VertexQ[CompleteGraph[4], 9]  (* 9 is not *)
```

```mathematica
In[1]:= VertexQ[CycleGraph[3], 1.0]  (* structural equality: 1.0 is not the vertex 1 *)
```

```mathematica
In[1]:= VertexQ[5, 2]  (* a non-graph first argument gives False, not unevaluated *)
```

### Notes

`VertexQ[g, v]` tests membership structurally, as `SameQ` does: the integer
vertex `1` is not matched by the real `1.0`. The test is an `O(1)` hash probe
into the graph's validated memo.

Unlike most graph accessors, `VertexQ` returns `False` (rather than staying
unevaluated) when the first argument is not a valid graph.

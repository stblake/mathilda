### Worked examples

```mathematica
In[1]:= GraphTriangleCount[CompleteGraph[4]]  (* C(4,3) = 4 triangles *)
```

```mathematica
In[1]:= GraphTriangleCount[CompleteGraph[5]]  (* C(5,3) = 10 triangles *)
```

```mathematica
In[1]:= GraphTriangleCount[CycleGraph[5]]  (* a cycle has none *)
```

```mathematica
In[1]:= GraphTriangleCount[PetersenGraph[]]  (* the Petersen graph is triangle-free *)
```

### Notes

`GraphTriangleCount[g]` counts triangles (undirected) or directed 3-cycles
`u->v->w->u` (directed). A complete graph `K_n` has `C(n, 3)` triangles.

The count is exact and uses a degree-ordered triangle listing that runs in
`O(m^1.5)`. `EdgeWeight` is ignored; a mixed graph is left unevaluated.

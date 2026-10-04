### Worked examples

```mathematica
In[1]:= EdgeList[GridGraph[{2, 3}]]  (* a 2 by 3 grid, first coordinate varying fastest *)
```

```mathematica
In[1]:= EdgeList[GridGraph[4]]  (* a bare integer gives a path on 4 vertices *)
```

```mathematica
In[1]:= VertexDegree[GridGraph[{3, 3}]]  (* corners 2, edge midpoints 3, centre 4 *)
```

```mathematica
In[1]:= VertexCount[GridGraph[{2, 3, 4}]]  (* a three-dimensional grid has the product of the dimensions *)
```

```mathematica
In[1]:= EdgeCount[GridGraph[{3, 4, 5}]]  (* the sum over axes of (n_i - 1) times the product of the others *)
```

```mathematica
In[1]:= EdgeList[GridGraph[{2, 2, 2}]]  (* the cube graph, identical to HypercubeGraph[3] *)
```

### Notes

`GridGraph[{n1, ..., nk}]` is the Cartesian product of paths of lengths `n1, ..., nk`. Vertex `1 + x1 + n1 x2 + n1 n2 x3 + ...` sits at coordinates `(x1, x2, ...)`, so the first coordinate runs fastest.

Dimensions must be positive integers (up to 32 of them). Very large grids, over 10^8 vertices or 5 x 10^7 edges, are refused and left unevaluated.

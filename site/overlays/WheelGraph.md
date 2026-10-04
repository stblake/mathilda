### Worked examples

```mathematica
In[1]:= EdgeList[WheelGraph[4]]  (* the smallest proper wheel is the complete graph on 4 vertices *)
```

```mathematica
In[1]:= EdgeList[WheelGraph[5]]  (* a hub joined to a 4-cycle rim *)
```

```mathematica
In[1]:= VertexDegree[WheelGraph[6]]  (* hub degree n - 1, every rim vertex degree 3 *)
```

```mathematica
In[1]:= EdgeCount[WheelGraph[10]]  (* 2 (n - 1) edges: spokes plus rim *)
```

```mathematica
In[1]:= EdgeList[WheelGraph[1]]  (* a single vertex and no edges *)
```

```mathematica
In[1]:= VertexCount[WheelGraph[7]]  (* the hub plus a rim of six vertices *)
```

### Notes

`WheelGraph[n]` has hub vertex 1 joined to every vertex of the cycle on `2..n`. For `n = 1` it is a single vertex; `n = 2` and `n = 3` would need repeated edges and are left unevaluated, as in Mathematica.

Vertex 1 always has degree `n - 1` and every rim vertex has degree 3.

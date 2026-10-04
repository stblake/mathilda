### Worked examples

```mathematica
In[1]:= EdgeList[KaryTree[7]]  (* a complete binary tree: 1 has children 2,3; 2 has 4,5; 3 has 6,7 *)
```

```mathematica
In[1]:= VertexCount[KaryTree[10]]  (* exactly n vertices *)
```

```mathematica
In[1]:= EdgeList[KaryTree[7, 3]]  (* ternary: vertex 1 has children 2,3,4 *)
```

```mathematica
In[1]:= VertexDegree[KaryTree[7]]  (* root and internal nodes have higher degree than the leaves *)
```

### Notes

`KaryTree[n]` is a binary tree on `n` vertices; `KaryTree[n, k]` a `k`-ary tree.
Vertices are numbered breadth-first, so vertex `i` has children `k(i-1)+2, ...,
k(i-1)+k+1` (1-based) up to `n`, and the last level fills left to right.

For the fully-filled tree given a number of *levels* rather than vertices, use
`CompleteKaryTree`.

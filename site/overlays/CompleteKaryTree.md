### Worked examples

```mathematica
In[1]:= EdgeList[CompleteKaryTree[3]]  (* three levels of a binary tree, 7 vertices *)
```

```mathematica
In[1]:= EdgeList[CompleteKaryTree[2, 3]]  (* a root with three children *)
```

```mathematica
In[1]:= VertexCount[CompleteKaryTree[4, 3]]  (* 1 + 3 + 9 + 27 vertices *)
```

```mathematica
In[1]:= VertexDegree[CompleteKaryTree[3, 2]]  (* root degree 2, internal nodes 3, leaves 1 *)
```

```mathematica
In[1]:= EdgeCount[CompleteKaryTree[5, 2]]  (* a tree has one fewer edge than vertices *)
```

```mathematica
In[1]:= GraphQ[CompleteKaryTree[3, 4]]  (* the result is an ordinary Graph expression *)
```

### Notes

The first argument counts **levels**, not vertices, so `CompleteKaryTree[n, k]` has `(k^n - 1)/(k - 1)` vertices. The second argument is the number of children per node and defaults to 2. Vertices are numbered in heap order: the children of vertex `i` are `k(i - 1) + 2 .. k(i - 1) + k + 1`.

Use `KaryTree[n, k]` instead when you want a given number of vertices with the last level partly filled.

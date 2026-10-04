### Worked examples

```mathematica
In[1]:= TreeGraphQ[StarGraph[5]]  (* a star is a tree *)
```

```mathematica
In[1]:= TreeGraphQ[CycleGraph[4]]  (* a cycle has as many edges as vertices *)
```

```mathematica
In[1]:= TreeGraphQ[PathGraph[{1, 2, 3, 4}]]  (* a path is a tree *)
```

```mathematica
In[1]:= TreeGraphQ[Graph[{1 -> 2, 1 -> 3}]]  (* direction is ignored, so an out-tree counts *)
```

```mathematica
In[1]:= TreeGraphQ[Graph[{1, 2, 3, 4}, {1 <-> 2, 3 <-> 4}]]  (* two components, so too few edges to be a tree *)
```

```mathematica
In[1]:= TreeGraphQ[Graph[{1 -> 2, 2 -> 1}]]  (* anti-parallel pair is two edges on two vertices *)
```

### Notes

`TreeGraphQ` ignores edge direction, so an out-tree counts, while the anti-parallel pair `1 -> 2, 2 -> 1` does not. The null graph is not a tree. The verdict is cached on the graph, so repeating the query on the same graph is `O(1)`.

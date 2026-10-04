### Worked examples

```mathematica
In[1]:= EdgeCount[CompleteGraph[5]]  (* K5 has n(n-1)/2 = 10 edges *)
```

```mathematica
In[1]:= VertexCount[CompleteGraph[4]]  (* vertices are 1..n *)
```

```mathematica
In[1]:= EdgeList[CompleteGraph[3]]  (* every pair joined, in row-major order *)
```

```mathematica
In[1]:= EdgeCount[CompleteGraph[{2, 3}]]  (* complete bipartite K(2,3): 2*3 = 6 edges *)
```

### Notes

`CompleteGraph[n]` joins every pair of the `n` vertices, so it has `n(n-1)/2`
edges and is the densest simple graph on `n` vertices.

`CompleteGraph[{n1, n2, ...}]` is the complete multipartite graph: the vertices
split into blocks of the given sizes and every pair in *different* blocks is
joined, none within a block. With two parts this is the complete bipartite
graph.

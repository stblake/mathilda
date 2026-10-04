### Worked examples

```mathematica
In[1]:= lg = LineGraph[PathGraph[4]];
```

```mathematica
In[1]:= {VertexList[lg], EdgeList[lg]}  (* three edges become three vertices along a path *)
```

### Notes

The line graph has one vertex per edge of the original graph, labelled by its `EdgeList`
position `1, 2, ..., m`. Two of these vertices are adjacent when the corresponding edges share an
endpoint (for a directed graph, when the head of one edge is the tail of the next).

So the line graph of a path `P_n` is a path `P_{n-1}`. A mixed graph (both directed and
undirected edges) is left unevaluated. The result is a canonical `Graph`; inspect it with
`VertexList` / `EdgeList`.

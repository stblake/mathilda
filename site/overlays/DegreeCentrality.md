### Worked examples

```mathematica
In[1]:= DegreeCentrality[StarGraph[5]]  (* the hub has degree 4, each leaf degree 1 *)
```

```mathematica
In[1]:= DegreeCentrality[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[1, 3], DirectedEdge[3, 2]}], "In"]  (* count incoming arcs only *)
```

```mathematica
In[1]:= DegreeCentrality[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[1, 3], DirectedEdge[3, 2]}], "Out"]  (* count outgoing arcs only *)
```

```mathematica
In[1]:= DegreeCentrality[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[1, 3], DirectedEdge[3, 2]}]]  (* the default is in-degree plus out-degree *)
```

```mathematica
In[1]:= DegreeCentrality[CompleteGraph[5]]  (* each vertex touches the other four *)
```

### Notes

The result is a list of exact integers in `VertexList` order. On an undirected graph every mode equals the ordinary vertex degree. On a mixed graph an undirected edge counts as one arc in each direction, so it adds 2 to the default total at each endpoint. Edge weights are ignored.

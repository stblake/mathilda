### Worked examples

```mathematica
In[1]:= EdgeList[EdgeAdd[PathGraph[{1, 2, 3}], 1 <-> 3]]  (* closing a path into a triangle *)
```

```mathematica
In[1]:= EdgeList[EdgeAdd[PathGraph[{1, 2, 3}], {3 <-> 4, 4 <-> 5}]]  (* a list of edges extends the path *)
```

```mathematica
In[1]:= VertexList[EdgeAdd[PathGraph[{1, 2, 3}], 3 <-> 9]]  (* an unseen endpoint becomes a new vertex *)
```

```mathematica
In[1]:= EdgeList[EdgeAdd[Graph[{1 -> 2, 2 -> 3}], 3 -> 1]]  (* in a directed graph the arrow keeps its direction *)
```

```mathematica
In[1]:= EdgeList[EdgeAdd[Graph[{1 <-> 2}], 2 -> 3]]  (* the arrow sugar is read as undirected here *)
```

```mathematica
In[1]:= EdgeList[EdgeAdd[Graph[{1, 2, 3}, {}], 1 <-> 2]]  (* edges can be added to an edgeless graph *)
```

### Notes

Graphs here are simple, so adding an edge that already exists (or a self-loop) is not accepted and `EdgeAdd` returns unevaluated. Endpoints that are not yet vertices are appended to the vertex list in order of first appearance.

Each edge may be written `u <-> v`, `u -> v`, `UndirectedEdge[u, v]` or `DirectedEdge[u, v]`. In a graph with no directed edges, `u -> v` is read as undirected.

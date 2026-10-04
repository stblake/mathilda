### Worked examples

```mathematica
In[1]:= EigenvectorCentrality[CycleGraph[5]]  (* every vertex of a cycle is equivalent *)
```

```mathematica
In[1]:= EigenvectorCentrality[StarGraph[5]]  (* the hub is twice as central as each leaf *)
```

```mathematica
In[1]:= EigenvectorCentrality[PathGraph[{1, 2, 3}]]  (* the middle vertex scores highest *)
```

```mathematica
In[1]:= EigenvectorCentrality[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3], DirectedEdge[3, 1], DirectedEdge[1, 3]}], "Out"]  (* scores from outgoing arcs *)
```

```mathematica
In[1]:= EigenvectorCentrality[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3], DirectedEdge[3, 1], DirectedEdge[1, 3]}], "In"]  (* scores from incoming arcs, the default *)
```

```mathematica
In[1]:= Total[EigenvectorCentrality[GridGraph[{3, 3}]]]  (* normalised to total one on a connected graph *)
```

### Notes

A vertex scores highly when it is linked to by other highly scoring vertices. The vector is the Perron eigenvector of the adjacency matrix, normalised to sum 1 on a connected graph. The optional `"In"` (default) or `"Out"` argument chooses whether incoming or outgoing arcs carry the score on a directed graph.

Disconnected and not strongly connected graphs are handled per strongly connected component, and vertices alone in a component score 0. The result is a list of machine reals.

### Worked examples

```mathematica
In[1]:= EdgeCoverQ[CycleGraph[4], {1 <-> 2, 3 <-> 4}]  (* two disjoint edges touch all four vertices *)
```

```mathematica
In[1]:= EdgeCoverQ[CycleGraph[4], {1 <-> 2}]  (* vertices 3 and 4 are missed *)
```

```mathematica
In[1]:= EdgeCoverQ[PathGraph[{1, 2, 3}], {2 <-> 1, 2 <-> 3}]  (* an undirected edge matches in either orientation *)
```

```mathematica
In[1]:= EdgeCoverQ[StarGraph[4], {UndirectedEdge[1, 2]}]  (* one spoke leaves the others uncovered *)
```

```mathematica
In[1]:= EdgeCoverQ[CycleGraph[4], {1 <-> 3}]  (* a pair that is not an edge of the graph gives false *)
```

```mathematica
In[1]:= EdgeCoverQ[x, {}]  (* a non-graph gives false *)
```

### Notes

The second argument is a list of edges of the graph, written `a <-> b` or `UndirectedEdge[a, b]` for undirected edges and `a -> b` for directed ones; a directed edge only matches as given. The predicate is `True` exactly when every vertex of the graph is an endpoint of one of the listed edges.

It is a linear check, never a search: use it to validate a candidate cover. A list element that is not an edge of the graph, or a non-list, gives `False`.

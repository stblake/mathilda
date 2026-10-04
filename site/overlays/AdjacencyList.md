### Worked examples

```mathematica
In[1]:= AdjacencyList[CycleGraph[4]]  (* each vertex of a 4-cycle has two neighbours *)
```

```mathematica
In[1]:= AdjacencyList[CompleteGraph[4], 1]  (* the neighbours of a single vertex *)
```

```mathematica
In[1]:= AdjacencyList[Graph[{1 -> 2, 1 -> 3, 2 -> 3}], 1]  (* directed: only successors count *)
```

```mathematica
In[1]:= AdjacencyList[Graph[{1 -> 2, 2 -> 3}], 3]  (* a directed sink has no successors *)
```

### Notes

Neighbours follow edge orientation: for a directed edge `v -> u` the successor `u` is a
neighbour of `v` (but `v` is not a neighbour of `u`), while an undirected edge makes each
endpoint a neighbour of the other. This matches the row convention of `AdjacencyMatrix`.

Neighbours are returned in first-appearance order and de-duplicated. `AdjacencyList[g]`
lists all vertices' neighbourhoods in `VertexList` order; `AdjacencyList[g, v]` gives just
one. A non-graph argument, or a `v` that is not a vertex of `g`, leaves the expression
unevaluated.

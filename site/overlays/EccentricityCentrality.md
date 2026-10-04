### Worked examples

```mathematica
In[1]:= EccentricityCentrality[CycleGraph[5]]  (* every vertex of a cycle has the same eccentricity *)
```

```mathematica
In[1]:= EccentricityCentrality[PathGraph[{1, 2, 3, 4, 5}]]  (* the middle vertex is the most central *)
```

```mathematica
In[1]:= EccentricityCentrality[StarGraph[5]]  (* the hub reaches every leaf in one step *)
```

```mathematica
In[1]:= EccentricityCentrality[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight -> {1, 2}]]  (* edge weights are used as path lengths *)
```

### Notes

The score of a vertex `v` is `1/VertexEccentricity[g, v]`, the reciprocal of the
distance to the farthest vertex `v` can reach. A vertex at the centre of the graph
has the smallest eccentricity and therefore the largest centrality; an isolated
vertex, whose eccentricity is `0`, is assigned the score `0`.

Distances come from a BFS on an unweighted graph and from a Dijkstra shortest path
when an `EdgeWeight` is supplied, so the centrality respects edge lengths. In a
disconnected graph the eccentricity is measured over the vertices actually
reachable from `v`, i.e. within its own component.

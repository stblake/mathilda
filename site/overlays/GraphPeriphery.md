### Worked examples

```mathematica
In[1]:= GraphPeriphery[PathGraph[{1, 2, 3, 4, 5}]]  (* the two endpoints are farthest out *)
```

```mathematica
In[1]:= GraphPeriphery[CycleGraph[5]]  (* a vertex-transitive graph: every vertex is peripheral *)
```

```mathematica
In[1]:= GraphPeriphery[CompleteGraph[4]]  (* all eccentricities equal 1 *)
```

### Notes

The periphery is the set of vertices whose eccentricity equals the graph
diameter — the vertices "on the rim". Dually, `GraphCenter` picks those of
minimum eccentricity (the radius).

On a vertex-transitive graph (a cycle, a complete graph) every vertex has the
same eccentricity, so the periphery is the whole vertex set. The periphery is
`{}` unless the graph is connected (strongly connected, if directed).

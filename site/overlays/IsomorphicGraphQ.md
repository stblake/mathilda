### Worked examples

```mathematica
In[1]:= IsomorphicGraphQ[CycleGraph[4], Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}]]  (* a relabelled 4-cycle *)
```

```mathematica
In[1]:= IsomorphicGraphQ[PathGraph[4], CycleGraph[4]]  (* different degree sequences *)
```

```mathematica
In[1]:= IsomorphicGraphQ[CompleteGraph[3], CycleGraph[3]]  (* the triangle is both *)
```

### Notes

`IsomorphicGraphQ[g1, g2, ...]` returns `True` when the graphs are all isomorphic — related by a
relabelling of vertices that preserves adjacency. It is variadic: all the arguments must be
mutually isomorphic.

Directed graphs are handled (in- and out-adjacency are matched separately), as are disconnected
graphs. A non-graph argument gives `False`.

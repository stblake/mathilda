### Worked examples

```mathematica
In[1]:= GraphDistance[CycleGraph[6], 1, 4]  (* opposite vertices of a 6-cycle are 3 hops apart *)
```

```mathematica
In[1]:= GraphDistance[PathGraph[{1, 2, 3, 4}], 1, 4]  (* a path's endpoints *)
```

```mathematica
In[1]:= GraphDistance[CycleGraph[5], 1]  (* two-argument form: all distances from a source *)
```

```mathematica
In[1]:= GraphDistance[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight -> {5, 7}], 1, 3]  (* weighted: Dijkstra, a machine real *)
```

```mathematica
In[1]:= GraphDistance[Graph[{1, 2, 3}, {1 <-> 2}], 1, 3]  (* an unreachable target is Infinity *)
```

### Notes

Without edge weights the distance is the BFS hop count, an integer; for a
directed graph the search follows edge direction. When the graph carries a
non-negative numeric `EdgeWeight`, the distance is the Dijkstra total and comes
back as a machine real (so `12.`, not `12`), agreeing bit-for-bit with the
single-source form and `GraphDistanceMatrix`.

`GraphDistance[g, s, t]` gives one length; `GraphDistance[g, s]` gives the list
of distances from `s` to every vertex in canonical order. An unreachable target
is `Infinity`; a symbolic or negative weight demotes the query to an unweighted
hop count rather than erroring.

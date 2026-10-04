### Worked examples

```mathematica
In[1]:= ClosenessCentrality[Graph[{1 <-> 2, 2 <-> 3}]]  (* the middle vertex is closest to the rest *)
```

```mathematica
In[1]:= ClosenessCentrality[Graph[{1 -> 2, 2 -> 3, 1 -> 3, 3 -> 4}]]  (* directed: only forward-reachable vertices count *)
```

### Notes

For each vertex `v` the value is `r/s`, where `r` is the number of vertices reachable from `v`
and `s` is the sum of their distances from `v`; a vertex that reaches nothing gets `0`.
Dividing by the reachable set (rather than by `n-1`) keeps the measure finite on a disconnected
or directed graph.

`EdgeWeight` is used as edge length when present, so on a weighted graph the distances — and
hence the centralities — reflect the weights rather than hop counts.

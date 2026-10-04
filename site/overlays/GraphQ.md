### Worked examples

```mathematica
In[1]:= GraphQ[CycleGraph[4]]  (* a generated graph is valid *)
```

```mathematica
In[1]:= GraphQ[5]  (* an atom is not a graph *)
```

```mathematica
In[1]:= GraphQ[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}]]  (* explicit vertices and edges *)
```

```mathematica
In[1]:= GraphQ[Graph[{1, 2}, {1 <-> 3}]]  (* an edge endpoint missing from the vertex list *)
```

```mathematica
In[1]:= GraphQ[{1 <-> 2, 2 <-> 3}]  (* a bare list of edges is not yet a graph *)
```

```mathematica
In[1]:= GraphQ[Graph[{1 -> 2, 2 -> 3}]]  (* directed graphs are valid too *)
```

### Notes

`GraphQ` is a structural test: it checks the `Graph[vertices, edges]` shape, that every edge endpoint is a listed vertex, and that edges are well formed. It never raises a message, and any non-graph gives `False`. The first call on a graph fills the validation memo that the other graph heads reuse.

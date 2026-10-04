### Worked examples

```mathematica
In[1]:= EdgeQ[CycleGraph[4], 1 <-> 2]  (* an undirected edge present *)
```

```mathematica
In[1]:= EdgeQ[CycleGraph[4], 2 <-> 1]  (* undirected edges match in either orientation *)
```

```mathematica
In[1]:= EdgeQ[CycleGraph[4], 1 <-> 3]  (* no such edge in the cycle *)
```

```mathematica
In[1]:= EdgeQ[Graph[{1 -> 2, 2 -> 3}], 1 -> 2]  (* a directed edge in its own orientation *)
```

```mathematica
In[1]:= EdgeQ[Graph[{1 -> 2, 2 -> 3}], 2 -> 1]  (* direction is never blurred *)
```

```mathematica
In[1]:= EdgeQ[Graph[{1 -> 2}], 1 <-> 2]  (* an undirected query does not match a directed edge *)
```

```mathematica
In[1]:= EdgeQ[5, 1 <-> 2]  (* a non-graph gives False *)
```

### Notes

`EdgeQ` accepts the edge sugar `u -> v` (directed) and `u <-> v` (undirected). An undirected query matches either orientation, but direction is never blurred between the two kinds. Membership is structural, so vertex `1` does not match `1.0`. Anything that is not a valid graph gives `False`.

### Worked examples

```mathematica
In[1]:= AcyclicGraphQ[Graph[{1 -> 2, 2 -> 3, 3 -> 4}]]  (* a directed chain is a DAG *)
```

```mathematica
In[1]:= AcyclicGraphQ[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]  (* the 1-2-3 cycle makes it False *)
```

```mathematica
In[1]:= AcyclicGraphQ[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4}]]  (* an undirected tree is a forest *)
```

```mathematica
In[1]:= AcyclicGraphQ[Graph[{1 -> 2, 2 -> 1}]]  (* anti-parallel edges are a 2-cycle *)
```

### Notes

A cycle is a closed walk using no edge twice, following directed edges forwards and undirected
edges either way. So an undirected graph is acyclic exactly when it is a forest, and a directed
graph exactly when it is a DAG. A pair `u -> v` together with `v -> u` counts as a cycle.

Being a predicate, it returns `False` (never stays unevaluated) for anything that is not an
acyclic graph, including a non-graph argument.

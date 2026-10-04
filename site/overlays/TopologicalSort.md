### Worked examples

```mathematica
In[1]:= TopologicalSort[Graph[{1 -> 2, 2 -> 3, 1 -> 3}]]  (* a DAG: u before v for every edge u -> v *)
```

```mathematica
In[1]:= TopologicalSort[{4 -> 2, 2 -> 1, 4 -> 3}]  (* the rule form builds Graph[...] first, then sorts *)
```

```mathematica
In[1]:= TopologicalSort[Graph[{1, 2, 3}, {}]]  (* an edgeless graph sorts to its VertexList *)
```

```mathematica
In[1]:= TopologicalSort[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]  (* a cycle has no topological order -- left unevaluated *)
```

```mathematica
In[1]:= TopologicalSort[CycleGraph[4]]  (* an undirected edge imposes no order, so this is not a DAG *)
```

### Notes

`TopologicalSort` is Kahn's algorithm: it repeatedly emits a vertex with no
remaining incoming edges. When several vertices are ready at once the tie is broken
by `VertexList` position (smallest index first), so the order is deterministic and
an already-sorted vertex list is returned unchanged.

The head applies only to a **directed acyclic graph**. A graph with a cycle, or one
carrying any undirected edge (an undirected edge imposes no `u` before `v`
ordering), is left unevaluated rather than given a best-effort order — so the
cyclic and undirected examples above print back as the input. A plain list of rules
is accepted as shorthand and built into a `Graph` with exactly the constructor's
edge sugar before sorting.

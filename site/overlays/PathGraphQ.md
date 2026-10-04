### Worked examples

```mathematica
In[1]:= PathGraphQ[PathGraph[{1, 2, 3, 4}]]  (* a simple path *)
```

```mathematica
In[1]:= PathGraphQ[CycleGraph[3]]  (* a cycle counts as a closed path *)
```

```mathematica
In[1]:= PathGraphQ[StarGraph[4]]  (* the hub has degree 3, so not a path *)
```

```mathematica
In[1]:= PathGraphQ[CompleteGraph[4]]  (* too many edges for every degree <= 2 *)
```

### Notes

`PathGraphQ[g]` is `True` when `g` is connected and every vertex has degree at
most 2 (in- and out-degree at most 1 for a directed graph). Following
Mathematica, a cycle satisfies this and so is a path; a vertex of degree 3 (as
in a star) or the density of a complete graph on 4+ vertices fails it.

A mixed graph is never a path.

### Worked examples

```mathematica
In[1]:= FindVertexCover[StarGraph[5]]  (* the centre alone covers every spoke *)
```

```mathematica
In[1]:= FindVertexCover[PathGraph[{a, b, c, d, e}]]  (* works on arbitrary vertex names *)
```

```mathematica
In[1]:= FindVertexCover[CompleteGraph[4]]  (* all but one vertex *)
```

```mathematica
In[1]:= FindVertexCover[CycleGraph[5]]  (* odd cycle needs three *)
```

```mathematica
In[1]:= VertexCoverQ[CycleGraph[5], FindVertexCover[CycleGraph[5]]]  (* the result passes the predicate *)
```

```mathematica
In[1]:= Length[FindVertexCover[PetersenGraph[]]]  (* minimum cover has size ten minus alpha *)
```

### Notes

The cover is a minimum one: it is the complement of a maximum independent set, found by the same exact solver as `FindIndependentVertexSet`. When several minimum covers exist the one returned is fixed by the vertex order, so repeated calls agree, but it need not be the one another system would pick.

Edge direction is ignored. Minimum vertex cover is NP-hard; if the search budget or a `TimeConstrained` limit is exhausted the call stays unevaluated rather than returning a non-minimum cover.

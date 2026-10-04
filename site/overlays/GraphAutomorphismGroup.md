### Worked examples

```mathematica
In[1]:= GraphAutomorphismGroup[PathGraph[{1, 2, 3, 4}]]  (* only the reversal *)
```

```mathematica
In[1]:= GraphAutomorphismGroup[CycleGraph[4]]  (* generators of the dihedral group of the square *)
```

```mathematica
In[1]:= GraphAutomorphismGroup[StarGraph[4]]  (* the leaves permute freely *)
```

```mathematica
In[1]:= Length[First[GraphAutomorphismGroup[PetersenGraph[]]]]  (* a handful of generators describes all 120 symmetries *)
```

```mathematica
In[1]:= GraphAutomorphismGroup[CompleteGraph[4]]  (* transpositions generate the symmetric group *)
```

### Notes

The group is returned as `PermutationGroup[{Cycles[...], ...}]` acting on the vertices numbered by their position (1-based) in the vertex list. Only a generating set is produced, and it is not unique, so do not compare groups by their generators.

A graph with no symmetry gives `PermutationGroup[{}]`. The engine is the individualization-refinement search also used by `FindGraphIsomorphism`, so large, sparse graphs are fast.

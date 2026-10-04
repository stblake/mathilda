### Worked examples

```mathematica
In[1]:= IndependentEdgeSetQ[CycleGraph[4], {1 <-> 2, 3 <-> 4}]  (* a matching: no shared vertex *)
```

```mathematica
In[1]:= IndependentEdgeSetQ[CycleGraph[4], {1 <-> 2, 2 <-> 3}]  (* both touch vertex 2 *)
```

### Notes

An independent edge set is a matching: a set of edges of `g` that pairwise share
no vertex. The test verifies every listed element is actually an edge of `g`
(an undirected edge matches either orientation) and that no vertex is used
twice.

A list containing something that is not an edge of `g` gives `False`. Use
`FindIndependentEdgeSet` to construct a maximum matching.

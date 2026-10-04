### Worked examples

```mathematica
In[1]:= EmptyGraphQ[Graph[{1, 2, 3}, {}]]  (* isolated vertices, no edges *)
```

```mathematica
In[1]:= EmptyGraphQ[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]  (* has edges *)
```

### Notes

A graph is empty when it has no edges, whatever its vertices — a set of isolated vertices is an
empty graph. The predicate is about edges only; `VertexCount` can be any value.

Like the other `...Q` predicates it returns `False` for a non-graph argument rather than staying
unevaluated.

### Worked examples

```mathematica
In[1]:= EdgeList[Subgraph[CompleteGraph[4], {1, 2, 3}]]  (* the induced triangle on three vertices *)
```

```mathematica
In[1]:= VertexCount[Subgraph[CycleGraph[5], {1, 2, 3}]]  (* keeps exactly the chosen vertices *)
```

```mathematica
In[1]:= EdgeList[Subgraph[CycleGraph[5], {1, 2, 3}]]  (* only edges with both endpoints kept *)
```

### Notes

`Subgraph[g, vs]` is the *vertex-induced* subgraph: it keeps the listed vertices
and exactly those edges of `g` whose endpoints both lie in the list. In a cycle,
inducing on three consecutive vertices keeps the two edges between them but drops
the wrap-around edge.

Vertices come out in the order given; list elements that are not vertices of `g`
are ignored. Edge weights are preserved.

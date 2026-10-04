### Worked examples

```mathematica
In[1]:= EdgeCount[CompleteGraph[4]]  (* a complete graph on 4 vertices has binomial(4,2) edges *)
```

```mathematica
In[1]:= EdgeCount[CycleGraph[5]]  (* a cycle has as many edges as vertices *)
```

```mathematica
In[1]:= EdgeCount[Graph[{1 -> 2, 2 -> 3}]]  (* directed edges count once each *)
```

### Notes

`EdgeCount[g]` is the length of `EdgeList[g]`: every edge of the canonical form,
directed or undirected, is counted exactly once.

The head takes a single argument and declines otherwise — there is no
edge-pattern counting form, so `EdgeCount[g, patt]` is left unevaluated. A
non-graph argument that is a `Hypergraph` is answered with its hyperedge count;
anything else stays unevaluated.

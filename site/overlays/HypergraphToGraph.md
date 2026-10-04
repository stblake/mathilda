### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= EdgeCount[HypergraphToGraph[h]]  (* directed edges v_a -> v_b for every a < b in a hyperedge *)
```

```mathematica
In[1]:= InputForm[HypergraphToGraph[{{1, 1, 2}}]]  (* a repeated vertex would be a self-loop, which is dropped *)
```

### Notes

`HypergraphToGraph` reads `h` as an **ordered** hypergraph (the Function
Repository convention): a hyperedge `{v1, ..., vk}` contributes the directed edge
`v_a -> v_b` for every position `a < b`. This is the one head besides the arity
family that cares about hyperedge order and repeats.

Because Mathilda graphs are simple, the result deviates from the FR multigraph: a
self-loop from a repeated vertex is dropped, parallel copies are merged, and every
vertex is kept — including one occurring only in a unary hyperedge, which the FR
function drops. The head also accepts a bare List of hyperedges.

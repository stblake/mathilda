### Worked examples

```mathematica
In[1]:= EdgeBetweennessCentrality[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4}]]  (* the central edge is on the most shortest paths *)
```

```mathematica
In[1]:= EdgeBetweennessCentrality[Graph[{1 -> 2, 2 -> 3, 1 -> 3, 3 -> 4}]]  (* in EdgeList order *)
```

### Notes

The result is a list of values in `EdgeList` order, one per edge: the sum over ordered pairs
of vertices `s`, `t` of the fraction of shortest `s`-`t` paths that use that edge.

Unlike the vertex `BetweennessCentrality`, this head uses `EdgeWeight` as edge lengths when the
graph carries usable weights, and it always sums over ordered pairs (there is no undirected ×1/2
factor). Tied shortest paths split the flow evenly.

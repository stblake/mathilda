### Worked examples

```mathematica
In[1]:= StronglyConnectedComponents[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]  (* the 3-cycle, then its sink *)
```

### Notes

Two vertices lie in the same strongly connected component when each can reach the other
following edge directions. The components come out in reverse-topological order of the
condensation (sinks first), the same order the directed `ConnectedComponents` uses.

An undirected edge is traversed in both directions, so on an undirected graph the strongly and
weakly connected components coincide.

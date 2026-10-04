### Worked examples

```mathematica
In[1]:= FindHamiltonianCycle[CycleGraph[5]]  (* the cycle itself, as a list of edges *)
```

```mathematica
In[1]:= FindHamiltonianCycle[PathGraph[4]]  (* a path has no Hamiltonian cycle *)
```

### Notes

The result is `{c}`, where `c` is a Hamiltonian cycle — a closed tour visiting every vertex
once — given as a list of edges in traversal order; it is `{}` when no such cycle exists.
`FindHamiltonianCycle[g, n]` returns up to `n` cycles and `[g, All]` returns all of them.

The search is exact backtracking with biconnectivity/strong-connectivity pruning, so both a
cycle and the answer `{}` are proofs. If the node budget is exhausted the call stays
unevaluated.

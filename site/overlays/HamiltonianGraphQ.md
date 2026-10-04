### Worked examples

```mathematica
In[1]:= HamiltonianGraphQ[CycleGraph[5]]  (* a cycle is its own Hamiltonian cycle *)
```

```mathematica
In[1]:= HamiltonianGraphQ[CompleteGraph[4]]  (* every complete graph on >= 3 vertices is Hamiltonian *)
```

```mathematica
In[1]:= HamiltonianGraphQ[PathGraph[{1, 2, 3}]]  (* a path has no cycle at all *)
```

```mathematica
In[1]:= HamiltonianGraphQ[PetersenGraph[]]  (* the classic non-Hamiltonian example *)
```

### Notes

`HamiltonianGraphQ[g]` tests for a Hamiltonian *cycle* (every vertex visited
once, returning to the start). The Petersen graph is the standard example of a
graph that is 3-regular and connected yet non-Hamiltonian.

The decision is exact, so `False` is a proof — but it is a backtracking search,
so a hard instance can exhaust the budget and leave the call unevaluated rather
than return a wrong answer; bound it with `TimeConstrained`.

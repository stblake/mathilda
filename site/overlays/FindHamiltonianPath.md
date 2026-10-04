### Worked examples

```mathematica
In[1]:= FindHamiltonianPath[CycleGraph[5]]  (* a cycle gives a path through every vertex *)
```

```mathematica
In[1]:= FindHamiltonianPath[PathGraph[{1, 2, 3}]]  (* the path itself *)
```

```mathematica
In[1]:= FindHamiltonianPath[CompleteGraph[4], 1, 3]  (* a Hamiltonian path from 1 to 3 *)
```

### Notes

The result is a list of vertices in visiting order, or `{}` when the graph has
no Hamiltonian path. `FindHamiltonianPath[g, s, t]` additionally fixes the two
endpoints.

The search is exact and complete, so `{}` is a proof of non-existence — but it
is backtracking, so a hard instance can exhaust the search budget and leave the
call unevaluated rather than return a wrong answer; bound it with
`TimeConstrained`.

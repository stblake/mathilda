### Worked examples

```mathematica
In[1]:= SeedRandom[42]; EdgeList[RandomHypergraph[{6, 3}, 2]]  (* 3 random 2-subsets of 1..6 *)
```

```mathematica
In[1]:= SeedRandom[42]; UniformHypergraphQ[RandomHypergraph[{6, 3}, 2], 2]  (* every hyperedge has arity 2, uniform by construction *)
```

```mathematica
In[1]:= SeedRandom[7]; HyperedgeSizes[RandomHypergraph[{10, {4, 3}}]]  (* the FR form: 4 hyperedges of arity 3 *)
```

### Notes

`RandomHypergraph[{n, m}, k]` draws `m` independent, uniformly random `k`-subsets
of `1..n` as hyperedges (so a hyperedge may repeat); the vertex list is all of
`1..n`, and the result is `k`-uniform by construction. `RandomHypergraph[{n, m}, k, c]`
gives `c` such hypergraphs.

`RandomHypergraph[{n, {e, a}}]` is the Function Repository form: `e` hyperedges of
arity `a` whose vertices are drawn **with replacement**, so a vertex may repeat
inside a hyperedge; `{n, {{e1, a1}, ...}}` mixes several arities. All forms draw
from the user-visible random stream, so `SeedRandom` reproduces them exactly — as
above, where the same seed makes the first two examples agree. Where the FR
function returns the bare edge List, this returns a validated `Hypergraph`;
`EdgeList` recovers the List.

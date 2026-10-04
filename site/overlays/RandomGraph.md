### Worked examples

```mathematica
In[1]:= r = RandomGraph[{6, 8}];
```

```mathematica
In[1]:= {VertexCount[r], EdgeCount[r]}  (* G(n, m): exactly n vertices and m edges *)
```

```mathematica
In[1]:= Length[RandomGraph[{4, 3}, 5]]  (* a list of k independent graphs *)
```

### Notes

`RandomGraph[{n, m}]` draws a uniformly random simple undirected graph on `n` vertices with
exactly `m` edges (the Erdős–Rényi `G(n, m)` model). The vertex and edge counts are therefore
fixed by the arguments; the structure is what varies, and it is reproducible under
`SeedRandom`.

`RandomGraph[{n, m}, k]` returns a list of `k` independent such graphs. `m` may not exceed
`n(n-1)/2`, the number of possible edges.

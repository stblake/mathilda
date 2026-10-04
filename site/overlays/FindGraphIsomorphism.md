### Worked examples

```mathematica
In[1]:= FindGraphIsomorphism[CycleGraph[4], CycleGraph[4]]  (* one map, here the identity *)
```

```mathematica
In[1]:= FindGraphIsomorphism[PathGraph[{1, 2, 3}], PathGraph[{b, a, c}]]  (* relabels the vertices of the first graph *)
```

```mathematica
In[1]:= FindGraphIsomorphism[PathGraph[{1, 2, 3}], PathGraph[{a, b, c}], 2]  (* at most two maps *)
```

```mathematica
In[1]:= Length[FindGraphIsomorphism[CycleGraph[4], CycleGraph[4], All]]  (* the dihedral group of the square has eight elements *)
```

```mathematica
In[1]:= FindGraphIsomorphism[CycleGraph[4], PathGraph[{1, 2, 3, 4}]]  (* non-isomorphic graphs give the empty list *)
```

### Notes

The result is a list of associations sending each vertex of the first graph to a vertex of the second; `{}` means the graphs are not isomorphic. A third argument `n` or `All` asks for up to that many maps, and `All` returns the whole coset of the automorphism group, which can be very large for symmetric graphs.

Mixed graphs, self-loops and multigraphs are supported; edge weights and other properties are ignored. If the search budget or a `TimeConstrained` limit runs out the call stays unevaluated rather than guessing.

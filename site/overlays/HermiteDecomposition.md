### Worked examples

```mathematica
In[1]:= HermiteDecomposition[{{1, 2}, {3, 4}}]  (* gives {u, r} with u unimodular and u . m == r *)
```

```mathematica
In[1]:= HermiteDecomposition[{{2, 3, 4}, {5, 6, 7}}]  (* a non-square 2x3 integer matrix *)
```

```mathematica
In[1]:= HermiteDecomposition[{{2, 0}, {0, 3}}]  (* already in Hermite form, so u is the identity *)
```

```mathematica
In[1]:= Det[HermiteDecomposition[{{1, 2}, {3, 4}}][[1]]]  (* the transform u has determinant +-1 *)
```

### Notes

`HermiteDecomposition[m]` returns `{u, r}` where `u` is unimodular
(`Abs[Det[u]] == 1`), `r` is the row Hermite normal form, and `u . m == r`. The
form `r` is in echelon shape with positive pivots and every entry above a pivot
reduced into `[0, pivot)`.

The decomposition is defined only over the integers. A non-integer matrix is left
unevaluated with a message, as is a non-rectangular or empty argument. The same
integer HNF primitive drives Mathilda's exact linear Diophantine solving.

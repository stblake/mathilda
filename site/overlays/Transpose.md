---
status: Stable
references:
  - "R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013 — the matrix transpose and index permutations of tensors."
---
### Worked examples

```mathematica
In[1]:= Transpose[{{1, 2, 3}, {4, 5, 6}}]
Out[1]= {{1, 4}, {2, 5}, {3, 6}}
```

```mathematica
In[1]:= Transpose[{{1, 2}, {3, 4}}]
Out[1]= {{1, 3}, {2, 4}}
```

```mathematica
In[1]:= Transpose[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {1, 1}]
Out[1]= {1, 5, 9}
```

The Gram matrix `A^T . A` of a symbolic `2x3` design matrix is symmetric by
construction:

```mathematica
In[1]:= Transpose[{{a, b, c}, {d, e, f}}] . {{a, b, c}, {d, e, f}}
Out[1]= {{a^2 + d^2, a b + d e, a c + d f}, {a b + d e, b^2 + e^2, b c + e f}, {a c + d f, b c + e f, c^2 + f^2}}
```

For an antisymmetric matrix `M = -M^T`, the sum `M + Transpose[M]` vanishes:

```mathematica
In[1]:= m = {{0, 1, 2}, {-1, 0, 3}, {-2, -3, 0}}; m + Transpose[m]
Out[1]= {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
```

### Notes

With one argument `Transpose` swaps the first two levels of a list, turning a `2x3` matrix into a `3x2` one. The optional permutation spec generalises this to arbitrary index reorderings of a rectangular array. A repeated index in the spec — `{1, 1}` in the third example — extracts the corresponding diagonal, here the main diagonal `{1, 5, 9}` of the `3x3` matrix. The spec must be a permutation of `{1, ..., r}` where `r` is the depth of the list.

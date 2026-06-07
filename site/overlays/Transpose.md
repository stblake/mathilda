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

### Notes

With one argument `Transpose` swaps the first two levels of a list, turning a `2x3` matrix into a `3x2` one. The optional permutation spec generalises this to arbitrary index reorderings of a rectangular array. A repeated index in the spec — `{1, 1}` in the third example — extracts the corresponding diagonal, here the main diagonal `{1, 5, 9}` of the `3x3` matrix. The spec must be a permutation of `{1, ..., r}` where `r` is the depth of the list.

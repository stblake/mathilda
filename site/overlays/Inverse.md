---
status: Stable
references:
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — Gauss-Jordan elimination and matrix inversion."
---
### Worked examples

```mathematica
In[1]:= Inverse[{{2, 0}, {0, 4}}]
Out[1]= {{1/2, 0}, {0, 1/4}}
```

```mathematica
In[1]:= Inverse[{{1, 2}, {3, 4}}]
Out[1]= {{-2, 1}, {3/2, -1/2}}
```

```mathematica
In[1]:= Inverse[{{1, 1, 1}, {0, 1, 1}, {0, 0, 1}}]
Out[1]= {{1, -1, 0}, {0, 1, -1}, {0, 0, 1}}
```

```mathematica
In[1]:= Inverse[{{1, 2}, {2, 4}}]
Out[1]= Inverse[{{1, 2}, {2, 4}}]
```

### Notes

Exact integer inputs yield exact rational inverses via fraction-free Gauss-Jordan elimination on the augmented matrix `[m | I]`. A singular matrix (last example) emits the `Inverse::sing` diagnostic on stderr and returns the call unevaluated; a non-square or empty argument emits `Inverse::matsq`. The `Method` option selects among `"DivisionFreeRowReduction"` (the default), `"OneStepRowReduction"`, and `"CofactorExpansion"`; an unrecognised method emits `Inverse::method`.

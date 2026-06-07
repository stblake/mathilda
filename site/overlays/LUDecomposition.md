---
status: Stable
references:
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — LU factorisation with partial pivoting."
  - "L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — Gaussian elimination and LU factorisation."
---
### Worked examples

```mathematica
In[1]:= LUDecomposition[{{4, 3}, {6, 3}}]
Out[1]= {{{4, 3}, {3/2, -3/2}}, {1, 2}, 0}
```

```mathematica
In[1]:= LUDecomposition[{{2, 1}, {4, 1}}]
Out[1]= {{{2, 1}, {2, -1}}, {1, 2}, 0}
```

The combined factor reconstructs `m` via `L . U` (here `L = {{1, 0}, {3/2, 1}}`, `U = {{4, 3}, {0, -3/2}}`):

```mathematica
In[1]:= {{1, 0}, {3/2, 1}} . {{4, 3}, {0, -3/2}}
Out[1]= {{4, 3}, {6, 3}}
```

### Notes

`LUDecomposition` returns a list `{lu, p, c}`: `lu` is the combined Doolittle factor whose strictly-lower triangle is `L` (with an implicit unit diagonal) and whose upper triangle is `U`; `p` is the 1-indexed row-permutation vector (here `{1, 2}`, i.e. no swap was needed); and `c` is an `L`-infinity condition estimate that is `0` for exact or symbolic input. The relation is `m[[p]] == L . U`, as the manual reconstruction above confirms. The algorithm is Doolittle elimination with partial pivoting; exact integer inputs keep exact rational factors. A singular matrix emits `LUDecomposition::sing` and completes with a zero pivot on `U`'s diagonal; a non-square or empty matrix emits `LUDecomposition::matsq`.

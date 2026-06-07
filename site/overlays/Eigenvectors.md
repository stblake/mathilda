---
status: Stable
references:
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — eigenvectors and invariant subspaces."
  - "R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013 — eigenspaces and diagonalisation."
---
### Worked examples

```mathematica
In[1]:= Eigenvectors[{{2, 1}, {0, 3}}]
Out[1]= {{1, 1}, {1, 0}}
```

```mathematica
In[1]:= Eigenvectors[{{2, 0}, {0, 5}}]
Out[1]= {{0, 1}, {1, 0}}
```

```mathematica
In[1]:= Eigenvectors[{{2, 1}, {1, 2}}]
Out[1]= {{1, 1}, {-1, 1}}
```

### Notes

The eigenvectors are listed in the same order as the corresponding eigenvalues from `Eigenvalues`, i.e. by decreasing absolute value of the eigenvalue. For exact or symbolic matrices the vectors come from the null-space pipeline and are **not** normalised — they are returned in a convenient integer/rational form (the symmetric example `{{1, 1}, {-1, 1}}` shows the orthogonal but unnormalised pair). For an `n x n` matrix the result always has length `n`; if the matrix is defective the shortfall is padded with zero vectors. Approximate numerical matrices instead return unit-`Norm` eigenvectors.

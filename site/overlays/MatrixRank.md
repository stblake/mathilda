---
status: Stable
references:
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — rank and the row echelon form."
  - "R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013 — rank and linear independence."
---
### Worked examples

```mathematica
In[1]:= MatrixRank[{{1, 2}, {2, 4}}]
Out[1]= 1
```

```mathematica
In[1]:= MatrixRank[{{1, 2}, {3, 4}}]
Out[1]= 2
```

```mathematica
In[1]:= MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= 2
```

### Notes

`MatrixRank` returns the number of linearly independent rows, which equals the number of independent columns. The first example is rank 1 because its rows are proportional, and the `3x3` example is rank 2 because its three rows satisfy a linear relation. The default exact path routes through `RowReduce` and counts the non-zero rows; rectangular matrices are accepted. For inexact (Real / MPFR) input, or with an explicit `Tolerance`, a numerical partial-pivot elimination is used and entries with `|entry| <= t` are treated as zero.

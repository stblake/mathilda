---
status: Stable
references:
  - "L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — Gram-Schmidt orthogonalisation and the QR factorisation."
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — QR factorisation algorithms."
---
### Worked examples

```mathematica
In[1]:= QRDecomposition[{{3, 0}, {0, 4}}]
Out[1]= {{{1, 0}, {0, 1}}, {{3, 0}, {0, 4}}}
```

```mathematica
In[1]:= QRDecomposition[{{1, 1}, {0, 1}}]
Out[1]= {{{1, 0}, {0, 1}}, {{1, 1}, {0, 1}}}
```

Reconstructing the original matrix from `{q, r}` uses `ConjugateTranspose[q] . r`:

```mathematica
In[1]:= q = QRDecomposition[{{1, 1}, {0, 1}}][[1]]; r = QRDecomposition[{{1, 1}, {0, 1}}][[2]]; ConjugateTranspose[q] . r
Out[1]= {{1, 1}, {0, 1}}
```

### Notes

`QRDecomposition` returns the **thin** factorisation `{q, r}` where `q` is row-orthonormal and `r` is upper triangular, with the original matrix recovered as `ConjugateTranspose[q] . r` (note the conjugate-transpose convention: `q`'s *rows* are the orthonormal basis). When the columns are already orthogonal — as for a diagonal matrix or the unit-column examples above — `q` is the identity and `r` reproduces the input. The algorithm is Modified Gram-Schmidt applied through the evaluator, so exact integer inputs keep `Sqrt[...]` column norms in closed form while machine-precision Real inputs return Real factors. `Pivoting -> True` adds a permutation matrix `p` so that `m . p == ConjugateTranspose[q] . r`.

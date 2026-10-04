### Worked examples

```mathematica
In[1]:= JordanDecomposition[{{2, 1}, {0, 2}}]  (* a defective matrix: j carries a 1 on the superdiagonal *)
```

```mathematica
In[1]:= JordanDecomposition[{{1, 0}, {0, 3}}]  (* distinct eigenvalues, so j is diagonal *)
```

```mathematica
In[1]:= JordanDecomposition[{{2, 1, 0}, {0, 2, 0}, {0, 0, 3}}]  (* one 2x2 Jordan block and one 1x1 block *)
```

### Notes

`JordanDecomposition[m]` returns `{s, j}` where `j` is the Jordan canonical form
of the square matrix `m` and `s` is the similarity matrix, so that
`m == s . j . Inverse[s]`. A `1` on `j`'s superdiagonal marks a deficient
eigenspace; the columns of `s` are the generalized eigenvectors grouped into the
Jordan chains that match `j`'s blocks.

Exact and symbolic matrices are decomposed from the characteristic polynomial and
the nullity sequence of `(m - lambda I)^k`, so the output stays exact. A generic
numeric matrix has distinct eigenvalues and is diagonalizable, so it takes the
numeric eigenvector path (`j` diagonal) and scales to large sizes; a numerically
defective matrix falls back to an exact decomposition and is numericalised back.
A non-square or empty matrix is left unevaluated.

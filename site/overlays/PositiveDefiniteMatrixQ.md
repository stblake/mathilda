### Worked examples

```mathematica
In[1]:= PositiveDefiniteMatrixQ[{{2, 1}, {1, 2}}]
Out[1]= True
```

```mathematica
In[1]:= PositiveDefiniteMatrixQ[{{1, 2}, {2, 1}}]
Out[1]= False
```

```mathematica
In[1]:= PositiveDefiniteMatrixQ[{{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}}]
Out[1]= True
```

```mathematica
In[1]:= PositiveDefiniteMatrixQ[{{2, I}, {-I, 2}}]
Out[1]= True
```

### Notes

`PositiveDefiniteMatrixQ[m]` tests whether the Hermitian part of `m` has only
positive eigenvalues, equivalently whether `Re[Conjugate[x] . m . x] > 0` for
every nonzero vector `x`. It is decided by attempting a Cholesky factorisation
of the Hermitian part. The second matrix `{{1,2},{2,1}}` has eigenvalues `3`
and `-1`, so it is indefinite and the test returns `False`. The third example
is the classic `3x3` second-difference (discrete Laplacian) matrix, which is
positive definite. The last example is a genuinely Hermitian complex matrix
`{{2, I}, {-I, 2}}` (note `Conjugate[I] = -I` off the diagonal) with
eigenvalues `1` and `3`; the complex path dispatches to LAPACK's `zpotrf` when
available. The function returns `False` on non-numeric, non-square, ragged, or
empty inputs.

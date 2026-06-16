### Worked examples

```mathematica
In[1]:= NegativeDefiniteMatrixQ[{{-2, 0}, {0, 3}}]
Out[1]= False
```

```mathematica
In[1]:= NegativeDefiniteMatrixQ[{{-3, 1, 0}, {1, -3, 1}, {0, 1, -3}}]
Out[1]= True
```

```mathematica
In[1]:= NegativeDefiniteMatrixQ[{{-2, I}, {-I, -2}}]
Out[1]= True
```

```mathematica
In[1]:= NegativeDefiniteMatrixQ[-Table[1/(i + j - 1), {i, 3}, {j, 3}]]
Out[1]= True
```

### Notes

`NegativeDefiniteMatrixQ[m]` tests whether `Re[Conjugate[x] . m . x] < 0` for
every nonzero `x`, equivalently that `-m` is positive definite. The diagonal
example fails because of the positive entry `3`. The symmetric tridiagonal matrix
is negative definite (eigenvalues all negative). The complex Hermitian matrix
`{{-2, I}, {-I, -2}}` shows that the test uses the Hermitian part. The last case
is the negated `3x3` Hilbert matrix — notoriously ill-conditioned yet still
detected as negative definite. The check is performed via an attempted Cholesky
factorisation of `-(m + ConjugateTranspose[m])/2`, dispatched to LAPACK on
numeric input.

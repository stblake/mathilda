### Worked examples

```mathematica
In[1]:= SchurDecomposition[{{1, 2}, {3, 4}}]  (* gives {q, t} with q orthonormal and t upper-triangular *)
```

```mathematica
In[1]:= SchurDecomposition[{{2, 0}, {0, 3}}]  (* an already-diagonal matrix, computed numerically *)
```

```mathematica
In[1]:= SchurDecomposition[N[{{3, -2}, {4, -1}}], RealBlockDiagonalForm -> False]  (* complex t, eigenvalues 1 +- 2 I on the diagonal *)
```

### Notes

`SchurDecomposition[m]` returns `{q, t}` with `q` orthonormal (unitary) and `t`
block upper-triangular, so that `m == q . t . ConjugateTranspose[q]`. By default
`RealBlockDiagonalForm -> True` keeps `t` real, using a 2×2 block for each
complex-conjugate eigenvalue pair; `-> False` makes `t` complex upper-triangular
with the eigenvalues on its diagonal. `Pivoting -> True` additionally returns a
scaling/permutation matrix.

The decomposition is numerical: a symbolic matrix has no closed-form Schur form
and is left unevaluated. `SchurDecomposition[{m, a}]` gives the generalized (QZ)
decomposition `{q, s, p, t}`. A machine or packed/`NDArray` matrix is read
straight off its buffer by the LAPACK kernel; arbitrary-precision real inputs use
an MPFR backend.

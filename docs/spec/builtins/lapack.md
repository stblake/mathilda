# LAPACK context

The ``LAPACK` `` context exposes the numerical **LAPACK** driver routines
directly in the REPL, under their canonical netlib names (`` LAPACK`dgesv ``,
`` LAPACK`dgeqrf ``, `` LAPACK`dgeev ``, …). These are the same double-precision
kernels Mathilda links for its machine-precision linear algebra (Apple
Accelerate on macOS, a system LAPACK on Linux).

They complement the high-level `System`` operators (`LinearSolve`, `Inverse`,
`Eigenvalues`, `SingularValueDecomposition`, …): use those for symbolic or
exact work, reach for ``LAPACK` `` when you want a specific numerical driver
and its exact factored outputs.

## Interface conventions

- **Ergonomic signatures.** Dimensions and leading dimensions are inferred from
  the array shapes; you pass only the operand arrays.
- **Array arguments** may be `NDArray[...]` values or nested `List`s. A
  right-hand side may be a vector or a matrix; the result matches its shape.
- **Multi-output routines return a `List`** of their outputs, e.g.
  `` LAPACK`dgeqrf[A] `` → `{Q, R}`, `` LAPACK`dgesdd[A] `` → `{U, S, VT}`,
  `` LAPACK`dgeev[A] `` → `{values, vectors}`.
- **Results** are `NDArray[...]` for real data and nested `List`s of
  `Complex[...]` for complex data. Apply `Normal` to a real matrix output
  before feeding it to `Transpose`, `DiagonalMatrix`, `Dot`, etc.
- **Unevaluated on failure.** A non-square input where one is required, a
  singular system, a non-positive-definite Cholesky input, a non-numeric entry,
  or non-convergence leaves the call unevaluated.
- All ``LAPACK` `` routines carry the `Protected` attribute and are present only
  when Mathilda is built with LAPACK (`USE_LAPACK`).

---

## Linear systems

### LAPACK`dgesv
`` LAPACK`dgesv[A, B] `` solves the square real system `A.X == B` by LU
factorisation (LAPACK `dgesv`). `B` is a length-`n` vector or an `n*k` matrix;
the result matches its shape.

```mathematica
In[1]:= LAPACK`dgesv[{{3, 2}, {1, 2}}, {7, 5}]
Out[1]= NDArray[{1.0, 2.0}]
```

### LAPACK`dtrtrs
`` LAPACK`dtrtrs[A, B] `` solves the triangular system `A.X == B`, treating `A`
as its upper-triangular part with a non-unit diagonal (LAPACK `dtrtrs`).

### LAPACK`dgels
`` LAPACK`dgels[A, B] `` gives the least-squares (`m >= n`) or minimum-norm
(`m < n`) full-rank solution of `A.X == B` for a real `m*n` matrix `A`
(LAPACK `dgels`).

```mathematica
In[1]:= LAPACK`dgels[{{1, 0}, {0, 1}, {1, 1}}, {1, 2, 4}]
Out[1]= NDArray[{1.5, 2.5}]
```

Complex counterparts: `` LAPACK`zgesv ``, `` LAPACK`ztrtrs ``, `` LAPACK`zgels ``.

---

## Factorizations

### LAPACK`dgetrf
`` LAPACK`dgetrf[A] `` gives `{LU, pivots}` for the real `m*n` matrix `A`: the
combined unit-lower/upper `LU` factors and the 1-based row pivots
(LAPACK `dgetrf`).

```mathematica
In[1]:= LAPACK`dgetrf[{{4, 3}, {6, 3}}]
Out[1]= {NDArray[{{6.0, 3.0}, {0.666667, 1.0}}], {2, 2}}
```

### LAPACK`dgetri
`` LAPACK`dgetri[A] `` gives the inverse of the square real matrix `A`, factored
internally with `dgetrf` then inverted with LAPACK `dgetri`. It is left
unevaluated for a singular matrix.

```mathematica
In[1]:= LAPACK`dgetri[{{4, 3}, {6, 3}}]
Out[1]= NDArray[{{-0.5, 0.5}, {1.0, -0.666667}}]
```

### LAPACK`dgeqrf
`` LAPACK`dgeqrf[A] `` gives the economy QR factorisation `{Q, R}` of the real
`m*n` matrix `A`: `Q` is `m*k`, `R` is `k*n` upper-trapezoidal, `k = Min[m, n]`,
and `A == Q.R` (LAPACK `dgeqrf` + `dorgqr`).

```mathematica
In[1]:= LAPACK`dgeqrf[{{1, 2}, {3, 4}}]
Out[1]= {NDArray[{{-0.316228, -0.948683}, {-0.948683, 0.316228}}],
         NDArray[{{-3.16228, -4.42719}, {0.0, -0.632456}}]}
```

### LAPACK`dpotrf
`` LAPACK`dpotrf[A] `` gives the lower Cholesky factor `L`
(`A == L.Transpose[L]`) of a real symmetric positive-definite matrix `A`
(LAPACK `dpotrf`). Unevaluated when `A` is not positive definite.

```mathematica
In[1]:= LAPACK`dpotrf[{{4, 2}, {2, 3}}]
Out[1]= NDArray[{{2.0, 0.0}, {1.0, 1.41421}}]
```

Complex counterparts: `` LAPACK`zgetrf ``, `` LAPACK`zgetri ``,
`` LAPACK`zgeqrf ``, `` LAPACK`zpotrf `` (with `A == L.ConjugateTranspose[L]`).

---

## Singular value decomposition

### LAPACK`dgesdd
`` LAPACK`dgesdd[A] `` gives the full SVD `{U, S, VT}` of a real `m*n` matrix
`A` via the divide-and-conquer driver (LAPACK `dgesdd`); `U` is `m*m`, `VT` is
`n*n`, and `S` has length `Min[m, n]`, with `A == U.DiagonalMatrix[S].VT` after
padding `S` to the right shape.

```mathematica
In[1]:= LAPACK`dgesdd[{{1, 2}, {3, 4}}]
Out[1]= {NDArray[{{-0.404554, -0.914514}, {-0.914514, 0.404554}}],
         NDArray[{5.46499, 0.365966}],
         NDArray[{{-0.576048, -0.817416}, {0.817416, -0.576048}}]}
```

### LAPACK`dgesvd
`` LAPACK`dgesvd[A] `` gives the full SVD `{U, S, VT}` via the QR-iteration
driver (LAPACK `dgesvd`).

Complex counterparts: `` LAPACK`zgesdd ``, `` LAPACK`zgesvd ``.

---

## Eigenproblems

### LAPACK`dgeev
`` LAPACK`dgeev[A] `` gives `{values, vectors}` for a general real square matrix
`A`: the eigenvalues (real or complex) and the matching right eigenvectors, so
that `A.vectors[[i]] == values[[i]] vectors[[i]]` (LAPACK `dgeev`).

```mathematica
In[1]:= LAPACK`dgeev[{{0, -1}, {1, 0}}]
Out[1]= {{0.0 + 1.0*I, 0.0 - 1.0*I},
         {{0.707107, 0.0 - 0.707107*I}, {0.707107, 0.0 + 0.707107*I}}}
```

### LAPACK`dsyev
`` LAPACK`dsyev[A] `` gives `{values, vectors}` for a real symmetric matrix `A`:
ascending real eigenvalues and orthonormal eigenvectors (upper triangle
referenced) (LAPACK `dsyev`).

```mathematica
In[1]:= LAPACK`dsyev[{{2, 1}, {1, 2}}]
Out[1]= {{1.0, 3.0}, {NDArray[{-0.707107, 0.707107}], NDArray[{0.707107, 0.707107}]}}
```

### LAPACK`zgeev / LAPACK`zheev
`` LAPACK`zgeev[A] `` is the complex general driver; `` LAPACK`zheev[A] `` is the
complex Hermitian driver (real ascending eigenvalues, unitary eigenvectors).

---

## Matrix norm

### LAPACK`dlange
`` LAPACK`dlange[A] `` gives the Frobenius norm of the real matrix `A`
(LAPACK `dlange`). `` LAPACK`zlange `` is the complex counterpart.

```mathematica
In[1]:= LAPACK`dlange[{{3, 4}}]
Out[1]= 5.0
```

---

See also the [LAPACK tutorial](../../../site/docs/tutorials/14-lapack.md) and the
[`BLAS` context](blas.md).

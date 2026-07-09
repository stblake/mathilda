# LAPACK drivers

This tutorial tours the ``LAPACK` `` context — direct access to the
machine-precision **LAPACK** driver routines for solving systems, factoring
matrices, computing the SVD, and solving eigenproblems. It builds on the
[BLAS tutorial](13-blas.md); the same conventions apply (arrays as lists or
`NDArray`, dimensions inferred, real results as `NDArray`).

Every transcript below was produced by the actual Mathilda binary.

## Solving linear systems

`dgesv` solves the square system `A.X == B` by LU factorisation. The
right-hand side can be a vector or a matrix; the result matches its shape.

```mathematica
In[1]:= LAPACK`dgesv[{{3, 2}, {1, 2}}, {7, 5}]
Out[1]= NDArray[{1.0, 2.0}]
```

For a rectangular `A`, `dgels` gives the least-squares (tall `A`) or
minimum-norm (wide `A`) solution:

```mathematica
In[2]:= LAPACK`dgels[{{1, 0}, {0, 1}, {1, 1}}, {1, 2, 4}]
Out[2]= NDArray[{1.5, 2.5}]
```

That is the best fit of `x` to the three equations `x1 = 1`, `x2 = 2`,
`x1 + x2 = 4` in the least-squares sense.

## Factorizations

`dgetrf` returns the LU factorisation as `{LU, pivots}` — the combined
unit-lower/upper factors and the 1-based row pivots:

```mathematica
In[3]:= LAPACK`dgetrf[{{4, 3}, {6, 3}}]
Out[3]= {NDArray[{{6.0, 3.0}, {0.666667, 1.0}}], {2, 2}}
```

`dgetri` inverts a matrix (factoring internally):

```mathematica
In[4]:= LAPACK`dgetri[{{4, 3}, {6, 3}}]
Out[4]= NDArray[{{-0.5, 0.5}, {1.0, -0.666667}}]
```

`dgeqrf` gives the economy QR factorisation `{Q, R}` with `A == Q.R`. Apply
`Normal` to the factors before recombining them, then `Chop` to drop the tiny
rounding residue:

```mathematica
In[5]:= Module[{q, r},
          {q, r} = Normal /@ LAPACK`dgeqrf[{{1, 2}, {3, 4}}];
          Chop[q . r]]
Out[5]= {{1.0, 2.0}, {3.0, 4.0}}
```

`dpotrf` gives the lower Cholesky factor `L` of a symmetric
positive-definite matrix, with `A == L.Transpose[L]`:

```mathematica
In[6]:= LAPACK`dpotrf[{{4, 2}, {2, 3}}]
Out[6]= NDArray[{{2.0, 0.0}, {1.0, 1.41421}}]
```

## Singular value decomposition

`dgesdd` (divide-and-conquer) and `dgesvd` (QR iteration) both return
`{U, S, VT}` with `A == U.DiagonalMatrix[S].VT`:

```mathematica
In[7]:= Module[{u, d, vt},
          {u, d, vt} = Normal /@ LAPACK`dgesdd[{{1, 2}, {3, 4}}];
          Chop[u . DiagonalMatrix[d] . vt]]
Out[7]= {{1.0, 2.0}, {3.0, 4.0}}
```

The singular values `S` come out in descending order.

## Eigenproblems

`dsyev` handles the real *symmetric* case, returning ascending eigenvalues and
orthonormal eigenvectors:

```mathematica
In[8]:= LAPACK`dsyev[{{2, 1}, {1, 2}}]
Out[8]= {{1.0, 3.0}, {NDArray[{-0.707107, 0.707107}], NDArray[{0.707107, 0.707107}]}}
```

`dgeev` handles a *general* matrix; its eigenvalues (and eigenvectors) may be
complex. A rotation matrix has purely imaginary eigenvalues:

```mathematica
In[9]:= LAPACK`dgeev[{{0, -1}, {1, 0}}]
Out[9]= {{0.0 + 1.0*I, 0.0 - 1.0*I},
         {{0.707107, 0.0 - 0.707107*I}, {0.707107, 0.0 + 0.707107*I}}}
```

The result is `{values, vectors}`, with `vectors[[i]]` the right eigenvector for
`values[[i]]` (so `A.vectors[[i]] == values[[i]] vectors[[i]]`).

For a complex Hermitian matrix, `zheev` returns real ascending eigenvalues:

```mathematica
In[10]:= LAPACK`zheev[{{2, Complex[0, 1]}, {Complex[0, -1], 2}}]
Out[10]= {{1.0, 3.0},
          {{-0.0 + 0.707107*I, -0.707107}, {0.0 + 0.707107*I, 0.707107}}}
```

## Matrix norm

`dlange` gives the Frobenius norm:

```mathematica
In[11]:= LAPACK`dlange[{{1, 2}, {3, 4}}]
Out[11]= 5.47723
```

## Complex counterparts and failures

Every real driver has a `z*` complex counterpart (`zgesv`, `zgetri`, `zgeqrf`,
`zgesdd`, `zgeev`, …). A singular system, a non-positive-definite Cholesky
input, or a symbolic entry leaves the call unevaluated rather than returning a
wrong answer:

```mathematica
In[12]:= LAPACK`dgetri[{{1, 2}, {2, 4}}]
Out[12]= LAPACK`dgetri[{{1, 2}, {2, 4}}]
```

(The matrix `{{1,2},{2,4}}` is singular.)

## Where to go next

The full reference is in
[`builtins/lapack.md`](../../../docs/spec/builtins/lapack.md); for the
lower-level BLAS kernels see the [BLAS tutorial](13-blas.md).

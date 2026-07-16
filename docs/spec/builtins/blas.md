# BLAS context

The ``BLAS` `` context exposes the numerical **BLAS** (Basic Linear Algebra
Subprograms) kernels directly in the REPL, under their canonical netlib names
(`` BLAS`ddot ``, `` BLAS`dgemv ``, `` BLAS`dgemm ``, …). These are the same
double-precision kernels Mathilda already links for its machine-precision
linear algebra (Apple Accelerate on macOS, a system BLAS on Linux); the
``BLAS` `` context lets you call the exact routine you want, with no symbolic
dispatch in the way.

They complement the high-level `System`` operators (`Dot`, `Norm`, …): use
those for symbolic or exact work, reach for ``BLAS` `` when you want a specific
machine-precision kernel — for benchmarking, teaching, or matching a numerical
recipe expressed in BLAS terms.

## Interface conventions

- **Ergonomic signatures.** Dimensions and leading dimensions are inferred from
  the array shapes, so you pass only the mathematically meaningful arguments —
  scalars (`alpha`, `beta`) and the operand arrays. There are no `m`, `n`, `k`,
  or `lda` arguments.
- **Array arguments** may be `NDArray[...]` values or ordinary nested `List`s
  (a vector is a flat list, a matrix a list of rows). Entries may be any exact
  or inexact real (for `d*` routines) or complex (for `z*` routines) number.
- **Results** come back as `NDArray[...]` for real routines and as a nested
  `List` of `Complex[...]` for complex routines (NDArray is real-only). Use
  `Normal` to turn a real result into an ordinary list.
- **Row-major.** The BLAS calls run in CBLAS row-major mode, matching
  `NDArray`'s native layout, so no transposition is needed.
- **Unevaluated on mismatch.** A shape mismatch, a symbolic entry, or (for a
  `d*` routine) a genuinely complex entry leaves the call unevaluated, so the
  expression flows through the evaluator unchanged.
- All ``BLAS` `` routines carry the `Protected` attribute. They are present
  only when Mathilda is built with BLAS/LAPACK (`USE_LAPACK`); otherwise the
  symbols are absent and the calls stay unevaluated.

---

## Level 1 — vector operations

### BLAS`ddot
`` BLAS`ddot[x, y] `` gives the dot product `x.y` of two equal-length real
vectors (BLAS `ddot`).

```mathematica
In[1]:= BLAS`ddot[{1, 2, 3}, {4, 5, 6}]
Out[1]= 32.0
```

### BLAS`dnrm2
`` BLAS`dnrm2[x] `` gives the Euclidean (2-)norm of the real vector `x`
(BLAS `dnrm2`).

```mathematica
In[1]:= BLAS`dnrm2[{3, 4}]
Out[1]= 5.0
```

### BLAS`dasum
`` BLAS`dasum[x] `` gives the sum of absolute values of the real vector `x`
(BLAS `dasum`).

### BLAS`idamax
`` BLAS`idamax[x] `` gives the **1-based** index of the entry of largest
absolute value in the real vector `x` (BLAS `idamax`).

### BLAS`daxpy
`` BLAS`daxpy[alpha, x, y] `` gives `alpha x + y` for a real scalar `alpha` and
equal-length real vectors `x`, `y` (BLAS `daxpy`).

```mathematica
In[1]:= BLAS`daxpy[2, {1, 1, 1}, {3, 4, 5}]
Out[1]= NDArray[{5.0, 6.0, 7.0}]
```

### BLAS`dscal
`` BLAS`dscal[alpha, x] `` gives `alpha x` for a real scalar `alpha` and real
vector `x` (BLAS `dscal`).

### Complex Level 1

- `` BLAS`zdotu[x, y] `` — unconjugated complex dot product `x.y` (`zdotu`).
- `` BLAS`zdotc[x, y] `` — conjugated complex dot product `Conjugate[x].y`
  (`zdotc`).
- `` BLAS`dznrm2[x] `` — Euclidean norm of a complex vector (`dznrm2`).
- `` BLAS`zaxpy[alpha, x, y] `` — `alpha x + y` with complex `alpha` (`zaxpy`).
- `` BLAS`zscal[alpha, x] `` — `alpha x` with complex `alpha` (`zscal`).

```mathematica
In[1]:= BLAS`zdotc[{Complex[1, 1]}, {Complex[1, 1]}]
Out[1]= 2.0

In[2]:= BLAS`zdotu[{I}, {I}]
Out[2]= -1.0
```

---

## Level 2 — matrix/vector operations

### BLAS`dgemv
`` BLAS`dgemv[alpha, A, x, beta, y] `` gives `alpha A.x + beta y` for a real
`m*n` matrix `A`, length-`n` vector `x`, and length-`m` vector `y`
(BLAS `dgemv`).

```mathematica
In[1]:= BLAS`dgemv[1, {{1, 2}, {3, 4}}, {1, 1}, 0, {0, 0}]
Out[1]= NDArray[{3.0, 7.0}]
```

### BLAS`dger
`` BLAS`dger[alpha, x, y, A] `` gives the rank-1 update
`alpha x.Transpose[{y}] + A` for a real scalar `alpha`, length-`m` vector `x`,
length-`n` vector `y`, and `m*n` matrix `A` (BLAS `dger`).

### BLAS`dtrmv
`` BLAS`dtrmv[A, x] `` gives `A.x`, treating `A` as the upper-triangular part of
the square real matrix `A` with a non-unit diagonal (BLAS `dtrmv`).

### BLAS`zgemv
`` BLAS`zgemv[alpha, A, x, beta, y] `` is the complex counterpart of
`` BLAS`dgemv `` (BLAS `zgemv`).

---

## Level 3 — matrix/matrix operations

### BLAS`dgemm
`` BLAS`dgemm[alpha, A, B, beta, C] `` gives `alpha A.B + beta C` for real
matrices with `A` of shape `m*k`, `B` of shape `k*n`, and `C` of shape `m*n`
(BLAS `dgemm`).

```mathematica
In[1]:= BLAS`dgemm[1, {{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}, 0, {{0, 0}, {0, 0}}]
Out[1]= NDArray[{{19.0, 22.0}, {43.0, 50.0}}]
```

With `beta` and a nonzero `C` it accumulates:

```mathematica
In[2]:= BLAS`dgemm[2, {{1, 0}, {0, 1}}, {{1, 2}, {3, 4}}, 1, {{1, 1}, {1, 1}}]
Out[2]= NDArray[{{3.0, 5.0}, {7.0, 9.0}}]
```

### BLAS`dsymm
`` BLAS`dsymm[alpha, A, B, beta, C] `` gives `alpha A.B + beta C` where `A` is
the `m*m` symmetric matrix given by its upper triangle, and `B`, `C` are `m*n`
(BLAS `dsymm`, left side).

### BLAS`dtrsm
`` BLAS`dtrsm[alpha, A, B] `` solves `A.X = alpha B` for `X`, where `A` is the
`m*m` upper-triangular (non-unit) matrix and `B` is `m*n` (BLAS `dtrsm`, left
side).

### BLAS`dsyrk
`` BLAS`dsyrk[alpha, A, beta, C] `` gives the symmetric rank-k update
`alpha A.Transpose[A] + beta C`, where `A` is `n*k` and `C` is the `n*n`
symmetric matrix given by its upper triangle (BLAS `dsyrk`). The full symmetric
matrix is returned.

### BLAS`zgemm
`` BLAS`zgemm[alpha, A, B, beta, C] `` is the complex counterpart of
`` BLAS`dgemm `` (BLAS `zgemm`).

### BLAS`zherk
`` BLAS`zherk[alpha, A, beta, C] `` gives the Hermitian rank-k update
`alpha A.ConjugateTranspose[A] + beta C` for real scalars `alpha`, `beta`, a
complex `n*k` matrix `A`, and a Hermitian `n*n` matrix `C` given by its upper
triangle (BLAS `zherk`). The full Hermitian matrix is returned.

---

## Unevaluated example

A symbolic entry (or a shape mismatch) leaves the call untouched:

```mathematica
In[1]:= BLAS`dgemm[1, {{a, b}}, {{c}, {d}}, 0, {{0}}]
Out[1]= BLAS`dgemm[1, {{a, b}}, {{c}, {d}}, 0, {{0}}]
```

See also the [BLAS tutorial](../../../site/docs/tutorials/13-blas.md) and the
[`LAPACK` context](lapack.md).

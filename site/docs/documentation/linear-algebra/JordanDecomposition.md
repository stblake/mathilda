# JordanDecomposition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`JordanDecomposition[m]`**

gives the Jordan decomposition of a square matrix m as a list {s, j}, where j is the Jordan canonical form of m and s is the similarity matrix, so that m == s . j . Inverse\[s\].  j is block-diagonal in Jordan blocks (an eigenvalue on the diagonal and 1's on the superdiagonal); the columns of s are the generalized eigenvectors grouped into Jordan chains that match j's blocks.  A 1 on j's superdiagonal marks a deficient eigenspace.

<details>
<summary>Notes</summary>

JordanDecomposition works on every input family supported by the rest of the linear-algebra builtins: - exact integer / rational matrices (output stays exact) - free-symbolic matrices with closed-form eigenvalues - complex matrices - machine-precision Real matrices (numeric eigensolver; the generic distinct-spectrum case is diagonalizable, so j is diagonal) - arbitrary-precision MPFR matrices (output at the input precision) For exact / symbolic input the eigenvalues come from the characteristic polynomial and the block structure from the nullity sequence of (m - lambda I)^k.  For inexact input a diagonalizable matrix takes the numeric eigenvector path; a numerically defective matrix is rationalised, decomposed exactly, and numericalised back.  An exact matrix with an irrational, defective eigenvalue whose generalized eigenspace cannot be spanned exactly is left unevaluated rather than returned wrong. A non-square or empty matrix emits JordanDecomposition::matsq and the call is left unevaluated.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= JordanDecomposition[{{27, 48, 81}, {-6, 0, 0}, {1, 0, 3}}][[2]]
Out[1]= {{6, 0, 0}, {0, 12, 1}, {0, 0, 12}}

In[2]:= {s, j} = JordanDecomposition[{{-103, -191, -255}, {110, 190, 222}, {9, 9, 33}}]; j
Out[2]= {{24, 0, 0}, {0, 48, 1}, {0, 0, 48}}

In[3]:= {s, j} = JordanDecomposition[{{-1.2, 2.7, 3.8}, {4.2, 4.4, 5.3}, {3.5, 7.6, 6.8}}]; Diagonal[j]
Out[3]= {13.7715, -2.12527, -1.6462}
```

## Algorithm

jordandecomp.c -- JordanDecomposition[m].

JordanDecomposition[m] gives {s, j} where j is the Jordan canonical form of the square matrix m and s is the similarity matrix, so that

```text
    m == s . j . Inverse[s].
```

j is block-diagonal in Jordan blocks J_k(lambda) (lambda on the diagonal, 1 on the superdiagonal); the columns of s are the generalized eigenvectors grouped into Jordan chains and ordered to match j's blocks.

------------------------------------------------------------------------ Two code paths, one chain engine.

```text
  Exact / symbolic (jd_exact_core):  the characteristic polynomial
  (eigen_char_poly_faddeev) is solved for the distinct eigenvalues with
  their algebraic multiplicities.  For each lambda, N = m - lambda*I, and
  the nullity sequence nu_i = dim ker(N^i) drives a top-down chain-top
  selection / bottom-up chain construction.  Every structural step runs in
  the field of the matrix entries (exact rational, or free-symbolic) using
  eigen_null_space + a MatrixRank-based "extend a spanning set to a basis"
  primitive.  A generalized eigenspace that cannot be spanned (an
  irrational *defective* eigenvalue, where is_zero_poly cannot decide the
  RowReduce pivots) leaves the whole call unevaluated rather than returning
  a wrong answer.

  Numeric (machine / MPFR, real / complex):  a generic numeric matrix has
  distinct eigenvalues and is diagonalizable, so jd_numeric_fast reads the
  eigenvectors straight off the numeric eigensolver (the Eigenvectors head,
  LAPACK-style Direct kernel / MPFR twin) as the columns of s and puts the
  per-column Rayleigh eigenvalue on the diagonal of j.  This is what every
  "diagonal j" example produces, and it scales (a 100x100 random matrix
  never touches exact arithmetic).  A numerically *defective* matrix (its
  eigenvectors are rank-deficient) falls back to rationalize -> exact core
  -> numericalize, which recovers the genuine block structure.
```

This file reuses the internal eigen helpers (eigen_char_poly_faddeev, eigen_solve_poly, eigen_null_space, ...) declared in eigen_internal.h. JordanDecomposition is part of the eigen/linalg cluster, so pulling in that header here is deliberate.

```text
Memory contract: standard builtin ownership (SPEC.md §4).  This file does
```

NOT free `res`; it returns a fresh Expr* the evaluator owns, or NULL to leave the call unevaluated.

## Implementation notes

- `Protected`.
- Works on every input family:
  - exact integer / rational matrices (exact `s` and `j`);
  - free-symbolic matrices with closed-form eigenvalues (e.g. any
    2×2, or small matrices with a factorable characteristic
    polynomial);
  - complex matrices;
  - machine-precision Real matrices (the numeric eigensolver; a generic
    distinct spectrum is diagonalizable, so `j` is diagonal);
  - arbitrary-precision MPFR matrices (output at the input precision).
- Algorithm:
  - **Exact / symbolic.** The eigenvalues come from the characteristic
    polynomial (`eigen_char_poly_faddeev` → `Solve`); for each eigenvalue
    `λ` the Jordan block structure is read from the nullity sequence
    `dim ker (m - λ I)^k`, and a top-down chain-top selection (via a
    `MatrixRank`-based "extend a spanning set to a basis") builds the
    Jordan chains bottom-up. The generalized eigenvectors of `s` are any
    valid Jordan basis, not a canonical one.
  - **Numeric.** A matrix with a distinct spectrum is diagonalizable, so
    `s` is the numeric eigenvectors as columns and `j = DiagonalMatrix`
    of the eigenvalues. A **machine-real** matrix (a packed / `NDArray`
    argument is read straight off its float64 buffer, no delist) uses a
    LAPACK-`dgeev` kernel that solves once and builds `s`/`j` directly
    from the raw buffers — no intermediate boxed eigenvector list, no
    second solve (a 100×100 in ~9.6 ms). A complex-entry matrix (or
    `USE_LAPACK=0`) instead runs one `Eigenvectors` solve and recovers
    each eigenvalue from its eigenvector by an `O(n)` component ratio;
    arbitrary-precision input keeps the full-precision solves. A matrix
    with a repeated numeric eigenvalue (defective, or ambiguous at
    machine precision) is rationalised, decomposed by the exact core, and
    numericalised back, which recovers the genuine block structure.
- Result-fidelity notes: the exact printed `s` is not canonical (any
  Jordan basis satisfying `m == s . j . Inverse[s]` is valid); numeric
  block-splitting is not bit-reproducible from Mathematica (an ill-posed
  problem), so the numeric contract is the residual `m ≈ s . j . Inverse[s]`
  plus a block structure matching the numeric nullities.
- Limitation: an exact matrix with an *irrational, defective* eigenvalue
  whose generalized eigenspace cannot be spanned exactly (the
  `is_zero_poly` pivot test cannot decide the row reduction) is left
  unevaluated rather than returned wrong.
- Issues `JordanDecomposition::argx` for the wrong argument count, and
  `JordanDecomposition::matsq` for a non-square or empty matrix; both
  leave the call unevaluated.
- Not lowered by `Compile[]` (the result is a heterogeneous pair of
  matrices, like `QRDecomposition` / `SingularValueDecomposition`); a
  packed / `NDArray[…]` argument is accepted (delisted and re-evaluated).

**Attributes:** `Protected`.

## References

**See also:** [Solve](../../solutions-of-equations/Solve/), [MatrixRank](../../linear-algebra/MatrixRank/), [NDArray](../../linear-algebra/NDArray/), [Eigenvectors](../../linear-algebra/Eigenvectors/), [QRDecomposition](../../linear-algebra/QRDecomposition/), [SingularValueDecomposition](../../linear-algebra/SingularValueDecomposition/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_jordandecomp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_jordandecomp.c)

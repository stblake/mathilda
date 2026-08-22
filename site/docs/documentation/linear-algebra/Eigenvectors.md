# Eigenvectors

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Eigenvectors[m]`**

gives a list of the eigenvectors of the square matrix m.

**`Eigenvectors[{m, a}]`**

gives the generalized eigenvectors of m with respect to a.

**`Eigenvectors[m, k]`**

gives the first k eigenvectors.

**`Eigenvectors[m, UpTo[k]]`**

gives k eigenvectors, or as many as are available.

<details>
<summary>Notes</summary>

For an n x n matrix Eigenvectors always returns a list of length n. If a matrix is defective for some eigenvalue, the corresponding shortfall is padded with zero vectors. For approximate numerical matrices the eigenvectors are normalised to unit Norm; for exact or symbolic matrices the eigenvectors are not normalised. Options: Cubics    -\> False      (radicals for a general cubic; True forces them) Quartics  -\> False      (radicals for a general quartic; True forces them) Method    -\> Automatic  (numeric-matrix method dispatch) With Cubics/Quartics -\> False (the default) a general irreducible cubic or quartic characteristic polynomial is returned as held Root\[\] objects. Special always-solvable families -- binomials, quadratic-in-x^m, and biquadratic-after-depression quartics -- are always returned in radical form regardless of these options. Method values for approximate-numeric matrices mirror Eigenvalues: Automatic, "Direct", "Arnoldi", "Banded", and "FEAST".  Each method returns the eigenvectors associated with the eigenvalues it would compute.  See ?Eigenvalues for the per-method semantics and sub-option grammar.  Non-numeric matrices ignore Method and use the symbolic null-space pipeline. Implementation status: "Direct" yields orthonormal eigenvectors for real symmetric matrices at machine precision (Householder + symmetric QR with accumulated rotations), unit-norm eigenvectors for real non-symmetric matrices via Hessenberg + Francis double- shift QR with accumulated Q followed by Schur-form back- substitution (complex eigenvalues yield complex eigenvectors emitted as Complex\[re, im\] entries), unitary orthonormal complex eigenvectors for complex Hermitian matrices via complex Householder tridiagonalisation + diagonal-phase correction + symmetric QR with composed complex Q, and unit-norm complex eigenvectors for complex non-Hermitian matrices via real block embedding into a 2n x 2n general matrix followed by grouped complex Gram-Schmidt extraction.  Automatic routes here. Arbitrary-precision (MPFR) inputs run a parallel "Direct" kernel at the input's combined precision: real symmetric (step 2d-A), real non-symmetric (step 2d-B), complex Hermitian (step 2d-C), and complex non-Hermitian (step 2d-D) MPFR all yield eigenvectors carrying full input precision -- orthonormal for the Hermitian / symmetric paths, unit 2-norm for the general paths. "Arnoldi" (Phase 3, machine + MPFR) returns Ritz vectors V\_m y\_i where V\_m is the orthonormal Arnoldi basis and y\_i diagonalises the small m x m Hessenberg H\_m.  Ritz vectors are unit 2-norm; for ill-conditioned matrices or m close to the spectral diameter they may need refinement (single inverse iteration is sufficient in practice).  MPFR Arnoldi carries input precision through to all output components. "Banded" (Phase 4, machine + MPFR) returns orthonormal real eigenvectors for real symmetric banded inputs and unitary complex eigenvectors for complex Hermitian banded inputs.  The band-Givens reduction accumulates an orthogonal (resp. unitary) Q during the chase; the final Z from the symmetric tridiag QR is composed against Q exactly as in the Direct Hermitian path. Banded refuses (falls back to Direct) on non-Hermitian or fully dense matrices. "FEAST" (Phase 5, machine + MPFR) returns the eigenvectors whose eigenvalues lie in the user-supplied real Interval -\> {a, b} -- orthonormal for real symmetric input, unitary for complex Hermitian input.  Sub-option grammar mirrors Eigen- values: Method -\> {"FEAST", "Interval" -\> {a, b}, "ContourPoints" -\> Ne, "SubspaceSize" -\> m0, "MaxIterations" -\> k, "Tolerance" -\> t}.  Same fail-soft cascade as Eigenvalues -- non-Hermitian, missing / degenerate Interval, generalised problem, Cholesky / LU failure, or non-convergence all fall back to Direct with a one-shot stderr warning.

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Eigenvectors[{{2, 1, 0}, {0, 2, 0}, {0, 0, 1}}]
Out[1]= {{1, 0, 0}, {0, 0, 0}, {0, 0, 1}}

In[2]:= Eigenvectors[{{1, 0, 1}, {0, 1, 0}, {0, 0, 1}}]
Out[2]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 0}}

In[3]:= Norm /@ Eigenvectors[{{1., 2.}, {2., 1.}}]
Out[3]= {1.0, 1.0}

In[4]:= FreeQ[Eigenvectors[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}], {0, 0, 0}]
Out[4]= True
```

### Applications (5)

```mathematica
In[5]:= Eigenvectors[{{2, 1}, {0, 3}}]
Out[5]= {{1, 1}, {1, 0}}

In[6]:= Eigenvectors[{{2, 0}, {0, 5}}]
Out[6]= {{0, 1}, {1, 0}}

In[7]:= Eigenvectors[{{2, 1}, {1, 2}}]
Out[7]= {{1, 1}, {-1, 1}}

In[8]:= Eigenvectors[{{a, b}, {c, d}}]
Out[8]= {{-b/(1/2 a - 1/2 d - 1/2 Sqrt[(-a - d)^2 - 4 (-b c + a d)]), 1}, {-b/(1/2 a - 1/2 d + 1/2 Sqrt[(-a - d)^2 - 4 (-b c + a d)]), 1}}

In[9]:= Eigenvectors[{{2, 0, 0}, {1, 2, 0}, {0, 1, 3}}]
Out[9]= {{0, 0, 1}, {0, -1, 1}, {0, 0, 0}}
```

## Implementation notes

**Algorithm.** `builtin_eigenvectors` shares its eigenvalue computation and `Method` dispatch with `Eigenvalues` (see that page: Faddeev–Leverrier characteristic polynomial for the exact/symbolic path; Householder tridiagonalisation + Wilkinson-shift symmetric QR, Hessenberg + Francis QR, and Banded/Arnoldi/FEAST kernels for the inexact path, with `LAPACK-HOOK` sites for `dtrevc`-style vector recovery). Once the eigenvalues are known and sorted by descending `Abs`, eigenvectors are obtained per eigenvalue by **null-space computation**: equal eigenvalues are collapsed into runs, the residual matrix `m − λI` (or `m − λa` for the generalised pencil; `a` itself for `Infinity` eigenvalues) is substituted and row-reduced, and `eigen_null_space` returns up to `multiplicity` basis vectors. A defective eigenvalue whose geometric multiplicity is short is padded in place with zero vectors so the `i`-th eigenvector lines up positionally with the `i`-th eigenvalue.

For inexact input the matrix is first **rationalised** (`common_rationalize_input` at the minimum precision present) so the rank defect needed to expose the eigenvector is not destroyed by floating-point round-off; the null-space and normalisation are done in exact rational arithmetic, then the result is numericalised back (`common_numericalize_result`) and each vector normalised to unit `Norm`. A `{k}`/`-k`/`UpTo[k]` spec trims the result.

**Data structures.** `Expr` matrices/vectors throughout the symbolic path; the residual substitution uses `ReplaceAll` of the internal `λ` symbol, and `eigen_null_space` drives the exact Gauss-Jordan row reducer. Numerical kernels use dense row-major `double`/MPFR buffers (shared with `Eigenvalues`).

**Complexity / limits.** Dominated by the eigenvalue solve plus one null-space (row-reduction) per distinct eigenvalue. Defective matrices yield fewer independent vectors than the eigenvalue multiplicity, made explicit by zero-vector padding; the generalised case is restricted to small `n`.

- `Protected`.
- For each eigenvalue `lambda_i` (with multiplicity μ), Eigenvectors
  computes the null space of `m - lambda_i I` via `RowReduce` and emits
  up to `μ` basis vectors. When the matrix is defective for that
  eigenvalue, the shortfall is padded in-line with zero vectors so the
  `i`-th eigenvector still corresponds positionally to the `i`-th
  eigenvalue.
- The returned list always has length `n` for an `n×n` matrix.
- **Irrational algebraic eigenvalues** bypass `RowReduce` entirely. Its
  pivot test (`is_zero_poly`) is a polynomial-*identity* test that treats
  each distinct radical as an independent generator, so for a *casus
  irreducibilis* eigenvalue — three real roots expressible in radicals
  only through complex cube roots — it cannot prove a genuinely-zero
  entry zero, and reports `m - lambda I` as full rank. When the null
  space comes back **empty** (never a legal answer: every eigenvalue has
  an eigenvector), the eigenvector is instead read off a column of
  `adj(m - x I)` reduced modulo the minimal polynomial `q` of `lambda`.
  The adjugate is built from cofactor determinants — no division, no
  pivoting — so nothing must be decided zero while it is computed, and
  the surviving "is this column zero?" test runs on a univariate
  polynomial over `Q`, where it is exact and complete. `q` is the
  characteristic polynomial itself when that is irreducible (checked
  with `IrreduciblePolynomialQ`), else `MinimalPolynomial[lambda, x]`.
  A *repeated* irrational eigenvalue drops the rank to `n-2` or below,
  vanishing the whole adjugate; the routine then declines rather than
  guessing, and the zero-vector padding above applies.
- For approximate matrices the result is computed in the rationalised
  domain, then numericalized and normalised to unit `Norm`. Exact /
  symbolic matrices return un-normalised eigenvectors.
- Generalised case: vectors that fall in the shared null space of `m`
  and `a` are returned as zero vectors (matching Mathematica's
  `Eigenvectors::geinsl1` warning behaviour).
- Options: same as Eigenvalues, including the numerical `Method`
  dispatch (`Automatic`, `"Direct"`, `"Arnoldi"`, `"Banded"`,
  `"FEAST"`).  `Method -> "FEAST"` returns the eigenvectors whose
  eigenvalues lie in the supplied `"Interval"` — orthonormal for
  real symmetric input, unitary for complex Hermitian input — and
  shares the same fail-soft cascade as `Eigenvalues`.

**Attributes:** `Protected`.

## References

**See also:** [RowReduce](../../linear-algebra/RowReduce/), [IrreduciblePolynomialQ](../../algebra/IrreduciblePolynomialQ/), [Norm](../../linear-algebra/Norm/), [Eigenvalues](../../linear-algebra/Eigenvalues/)

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — eigenvectors and invariant subspaces.
- R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013 — eigenspaces and diagonalisation.
- Gene H. Golub, Charles F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins University Press, 2013).
- J. H. Wilkinson, *The Algebraic Eigenvalue Problem* (Oxford University Press, 1965).
- Source: [`src/linalg/eigen.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/eigen.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_eigen.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eigen.c)
- Tests: [`tests/test_mateigen_arnoldi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mateigen_arnoldi.c)
- Tests: [`tests/test_mateigen_banded.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mateigen_banded.c)
- Tests: [`tests/test_mateigen_direct.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mateigen_direct.c)

## Notes & additional examples

### Notes

The eigenvectors are listed in the same order as the corresponding eigenvalues from `Eigenvalues`, i.e. by decreasing absolute value of the eigenvalue. For exact or symbolic matrices the vectors come from the null-space pipeline and are **not** normalised — they are returned in a convenient integer/rational form (the symmetric example `{{1, 1}, {-1, 1}}` shows the orthogonal but unnormalised pair). For an `n x n` matrix the result always has length `n`; if the matrix is defective the shortfall is padded with zero vectors. Approximate numerical matrices instead return unit-`Norm` eigenvectors.

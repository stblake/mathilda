# Eigenvectors

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Eigenvectors[m]
    gives a list of the eigenvectors of the square matrix m.
Eigenvectors[{m, a}]
    gives the generalized eigenvectors of m with respect to a.
Eigenvectors[m, k]
    gives the first k eigenvectors.
Eigenvectors[m, UpTo[k]]
    gives k eigenvectors, or as many as are available.

For an n x n matrix Eigenvectors always returns a list of length n.
If a matrix is defective for some eigenvalue, the corresponding
shortfall is padded with zero vectors. For approximate numerical
matrices the eigenvectors are normalised to unit Norm; for exact
or symbolic matrices the eigenvectors are not normalised.

Options:
    Cubics    -> True       (use radicals to solve cubics)
    Quartics  -> True       (use radicals to solve quartics)
    Method    -> Automatic  (numeric-matrix method dispatch)

Method values for approximate-numeric matrices mirror Eigenvalues:
Automatic, "Direct", "Arnoldi", "Banded", and "FEAST".  Each
method returns the eigenvectors associated with the eigenvalues
it would compute.  See ?Eigenvalues for the per-method semantics
and sub-option grammar.  Non-numeric matrices ignore Method and
use the symbolic null-space pipeline.

Implementation status: "Direct" yields orthonormal eigenvectors
for real symmetric matrices at machine precision (Householder +
symmetric QR with accumulated rotations), unit-norm eigenvectors
for real non-symmetric matrices via Hessenberg + Francis double-
shift QR with accumulated Q followed by Schur-form back-
substitution (complex eigenvalues yield complex eigenvectors
emitted as Complex[re, im] entries), unitary orthonormal complex
eigenvectors for complex Hermitian matrices via complex
Householder tridiagonalisation + diagonal-phase correction +
symmetric QR with composed complex Q, and unit-norm complex
eigenvectors for complex non-Hermitian matrices via real block
embedding into a 2n x 2n general matrix followed by grouped
complex Gram-Schmidt extraction.  Automatic routes here.
Arbitrary-precision (MPFR) inputs run a parallel "Direct" kernel
at the input's combined precision: real symmetric (step 2d-A),
real non-symmetric (step 2d-B), complex Hermitian (step 2d-C),
and complex non-Hermitian (step 2d-D) MPFR all yield eigenvectors
carrying full input precision -- orthonormal for the Hermitian /
symmetric paths, unit 2-norm for the general paths.
"Arnoldi" (Phase 3, machine + MPFR) returns Ritz vectors
V_m y_i where V_m is the orthonormal Arnoldi basis and y_i
diagonalises the small m x m Hessenberg H_m.  Ritz vectors are
unit 2-norm; for ill-conditioned matrices or m close to the
spectral diameter they may need refinement (single inverse
iteration is sufficient in practice).  MPFR Arnoldi carries
input precision through to all output components.
"Banded" (Phase 4, machine + MPFR) returns orthonormal real
eigenvectors for real symmetric banded inputs and unitary
complex eigenvectors for complex Hermitian banded inputs.  The
band-Givens reduction accumulates an orthogonal (resp. unitary)
Q during the chase; the final Z from the symmetric tridiag QR
is composed against Q exactly as in the Direct Hermitian path.
Banded refuses (falls back to Direct) on non-Hermitian or fully
dense matrices.
"FEAST" (Phase 5, machine + MPFR) returns the eigenvectors
whose eigenvalues lie in the user-supplied real Interval ->
{a, b} -- orthonormal for real symmetric input, unitary for
complex Hermitian input.  Sub-option grammar mirrors Eigen-
values: Method -> {"FEAST", "Interval" -> {a, b},
"ContourPoints" -> Ne, "SubspaceSize" -> m0,
"MaxIterations" -> k, "Tolerance" -> t}.  Same fail-soft
cascade as Eigenvalues -- non-Hermitian, missing / degenerate
Interval, generalised problem, Cholesky / LU failure, or
non-convergence all fall back to Direct with a one-shot
stderr warning.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Eigenvectors[{{2, 1, 0}, {0, 2, 0}, {0, 0, 1}}]
Out[1]= {{1, 0, 0}, {0, 0, 0}, {0, 0, 1}}

In[2]:= Eigenvectors[{{1, 0, 1}, {0, 1, 0}, {0, 0, 1}}]
Out[2]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 0}}

In[3]:= Norm /@ Eigenvectors[{{1., 2.}, {2., 1.}}]
Out[3]= {1.0, 1.0}

In[4]:= Norm /@ Eigenvectors[{{1, 2}, {2, 1}}]
Out[4]= {Sqrt[2], Sqrt[2]}
```

## Implementation notes

**Algorithm.** `builtin_eigenvectors` shares its eigenvalue computation and `Method` dispatch with `Eigenvalues` (see that page: Faddeev–Leverrier characteristic polynomial for the exact/symbolic path; Householder tridiagonalisation + Wilkinson-shift symmetric QR, Hessenberg + Francis QR, and Banded/Arnoldi/FEAST kernels for the inexact path, with `LAPACK-HOOK` sites for `dtrevc`-style vector recovery). Once the eigenvalues are known and sorted by descending `Abs`, eigenvectors are obtained per eigenvalue by **null-space computation**: equal eigenvalues are collapsed into runs, the residual matrix `m − λI` (or `m − λa` for the generalised pencil; `a` itself for `Infinity` eigenvalues) is substituted and row-reduced, and `eigen_null_space` returns up to `multiplicity` basis vectors. A defective eigenvalue whose geometric multiplicity is short is padded in place with zero vectors so the `i`-th eigenvector lines up positionally with the `i`-th eigenvalue.

For inexact input the matrix is first **rationalised** (`common_rationalize_input` at the minimum precision present) so the rank defect needed to expose the eigenvector is not destroyed by floating-point round-off; the null-space and normalisation are done in exact rational arithmetic, then the result is numericalised back (`common_numericalize_result`) and each vector normalised to unit `Norm`. A `{k}`/`-k`/`UpTo[k]` spec trims the result.

**Data structures.** `Expr` matrices/vectors throughout the symbolic path; the residual substitution uses `ReplaceAll` of the internal `λ` symbol, and `eigen_null_space` drives the exact Gauss-Jordan row reducer. Numerical kernels use dense row-major `double`/MPFR buffers (shared with `Eigenvalues`).

**Complexity / limits.** Dominated by the eigenvalue solve plus one null-space (row-reduction) per distinct eigenvalue. Defective matrices yield fewer independent vectors than the eigenvalue multiplicity, made explicit by zero-vector padding; the generalised case is restricted to small `n`.

- `Protected`.
- For each eigenvalue `lambda_i` (with multiplicity μ), Eigenvectors

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — eigenvectors and invariant subspaces.
- R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013 — eigenspaces and diagonalisation.
- Gene H. Golub, Charles F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins University Press, 2013).
- J. H. Wilkinson, *The Algebraic Eigenvalue Problem* (Oxford University Press, 1965).
- Source: [`src/linalg/eigen.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/eigen.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Eigenvectors[{{2, 1}, {0, 3}}]
Out[1]= {{1, 1}, {1, 0}}
```

```mathematica
In[1]:= Eigenvectors[{{2, 0}, {0, 5}}]
Out[1]= {{0, 1}, {1, 0}}
```

```mathematica
In[1]:= Eigenvectors[{{2, 1}, {1, 2}}]
Out[1]= {{1, 1}, {-1, 1}}
```

```mathematica
In[1]:= Eigenvectors[{{a, b}, {c, d}}]
Out[1]= {{-b/(1/2 a - 1/2 d - 1/2 Sqrt[(-a - d)^2 - 4 (-b c + a d)]), 1}, {-b/(1/2 a - 1/2 d + 1/2 Sqrt[(-a - d)^2 - 4 (-b c + a d)]), 1}}
```

```mathematica
In[1]:= Eigenvectors[{{2, 0, 0}, {1, 2, 0}, {0, 1, 3}}]
Out[1]= {{0, 0, 1}, {0, -1, 1}, {0, 0, 0}}
```

### Notes

The eigenvectors are listed in the same order as the corresponding eigenvalues from `Eigenvalues`, i.e. by decreasing absolute value of the eigenvalue. For exact or symbolic matrices the vectors come from the null-space pipeline and are **not** normalised — they are returned in a convenient integer/rational form (the symmetric example `{{1, 1}, {-1, 1}}` shows the orthogonal but unnormalised pair). For an `n x n` matrix the result always has length `n`; if the matrix is defective the shortfall is padded with zero vectors. Approximate numerical matrices instead return unit-`Norm` eigenvectors.
